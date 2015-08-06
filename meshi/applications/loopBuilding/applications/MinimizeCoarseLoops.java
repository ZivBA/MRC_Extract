package meshi.applications.loopBuilding.applications;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.StringTokenizer;

import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.compositeTorsions.ramachandran.RamachandranCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBondEnergy;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_Minimization_BBonly_Creator;
import meshi.energy.softExcludedVol.SoftExcludedVolCreator;
import meshi.energy.tether.TetherCreator;
import meshi.energy.tether.TetherEnergy;
import meshi.energy.torsionVal.TorsionValCreator;
import meshi.energy.torsionVal.TorsionValEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LBFGS;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.Minimizer;
import meshi.optimizers.MinimizerException;
import meshi.optimizers.SteepestDecent;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;
import programs.PutHydrogens;


public class MinimizeCoarseLoops extends MeshiProgram implements Residues, AtomTypes{ 

	private static CommandList commands; 
	private static String commandsFileName = null;
	private static String modelFileName = null;  
	private static String loopsFileName = null;  	
	private static String dataFileName = null;  
	private static String noTorsionFileName = null;  
	private static String refFileName = null;  
	private static int resStart = -999;
	private static int resEnd = -999;


	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		Protein reference = null;
		Protein model = null;
		DistanceMatrix distanceMatrix = null;
		TotalEnergy energy = null;
		Minimizer minimizer = null;
		int dontTakeTorsion1 = -999;
		int dontTakeTorsion2 = -999;


		// Not taking the torsions:
		String[] closures = File2StringArray.f2a(noTorsionFileName);
		if (closures.length==2) {
			StringTokenizer st = new StringTokenizer(closures[0]);
			dontTakeTorsion1 = (new Integer(st.nextToken())).intValue(); 
			st = new StringTokenizer(closures[1]);
			dontTakeTorsion2 = (new Integer(st.nextToken())).intValue(); 
		}
		else 
			throw new RuntimeException("There should be only 2 closure residues in the file: " + noTorsionFileName);
		
		
		// The creators for the terms
		EnergyCreator[] energyCreators = {  
				new BondCreator(),
				new AngleCreator(),
				new PlaneCreator(),
				new OutOfPlaneCreator(),
				new SoftExcludedVolCreator(10.0 , 12 , 1.0),
				new SimpleHydrogenBond_Dahiyat_Minimization_BBonly_Creator(3.0),
				new TorsionValCreator(2.0),
				new TorsionValCreator(10*2.0),
				new RamachandranCreator(1.5),
				new TetherCreator(0.0000000001, new AtomList.BackboneFilter())
		};	

		// Loading the reference and the background model	
//		reference = new Protein((new AtomList(refFileName)).noOXTFilter().backbone(), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
//		model = new Protein((new AtomList(modelFileName)).noOXTFilter().backbone(), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		reference = new Protein((new AtomList(refFileName)).noOXTFilter().backbone(), new ResidueExtendedAtoms(ADD_HYDROGENS_AND_FREEZE));
		PutHydrogens.adjustHydrogens(commands, reference);
		model = new Protein((new AtomList(modelFileName)).noOXTFilter().backbone(), new ResidueExtendedAtoms(ADD_HYDROGENS_AND_FREEZE));
		PutHydrogens.adjustHydrogens(commands, model);
		for (int cc=0 ; cc<model.atoms().size() ; cc++)
			model.atoms().atomAt(cc).setChain("A");
		model.freeze();
		for (int c=resStart ; c<=resEnd ; c++)
			model.residue(c).atoms().defrost();
		distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
		energy = new TotalEnergy(model, distanceMatrix, energyCreators, commands);
		// Setting the tether term reliabilities
		for (int atomC=0 ; atomC<model.residue(dontTakeTorsion1).atoms().size() ; atomC++) {
			if (model.residue(dontTakeTorsion1).atoms().atomAt(atomC).reliability()>0.99) {
				model.residue(dontTakeTorsion1).atoms().atomAt(atomC).setReliability(0.25);
				System.out.println("Reliability: " + model.residue(dontTakeTorsion1).atoms().atomAt(atomC).residueNumber() + " " + 
						model.residue(dontTakeTorsion1).atoms().atomAt(atomC).name() + ": " + model.residue(dontTakeTorsion1).atoms().atomAt(atomC).reliability());
			}
		}
		for (int atomC=0 ; atomC<model.residue(dontTakeTorsion2).atoms().size() ; atomC++) {
			if (model.residue(dontTakeTorsion2).atoms().atomAt(atomC).reliability()>0.99) {
				model.residue(dontTakeTorsion2).atoms().atomAt(atomC).setReliability(0.25);
				System.out.println("Reliability: " + model.residue(dontTakeTorsion2).atoms().atomAt(atomC).residueNumber() + " " + 
						model.residue(dontTakeTorsion2).atoms().atomAt(atomC).name() + ": " + model.residue(dontTakeTorsion2).atoms().atomAt(atomC).reliability());
			}
		}
		/* This was an older setting which was not working well
		// Setting the tether term reliabilities
		double halfLoopLength = (resEnd-resStart+1)/2.0 - 1.0;
		for (int res=resStart ; res<=resEnd ; res++) {
			for (int atomC=0 ; atomC<model.residue(res).atoms().size() ; atomC++) {
				if (model.residue(res).atoms().atomAt(atomC).reliability()>0.99) {
					int disFromEnd = Math.min(res-resStart, resEnd-res);
					model.residue(res).atoms().atomAt(atomC).setReliability(Math.max(0.05 , 1.0-disFromEnd*0.95/halfLoopLength));
					System.out.println("Reliability: " + model.residue(res).atoms().atomAt(atomC).residueNumber() + " " + 
							model.residue(res).atoms().atomAt(atomC).name() + ": " + model.residue(res).atoms().atomAt(atomC).reliability());
				}
			}
		}
		*/
		SimpleHydrogenBondEnergy hbTerm = (SimpleHydrogenBondEnergy) energy.getEnergyTerm(new SimpleHydrogenBondEnergy());
		TorsionValEnergy torsionTerm1 = (TorsionValEnergy) energy.getEnergyTerms(new TorsionValEnergy())[0];
		TorsionValEnergy torsionTerm2 = (TorsionValEnergy) energy.getEnergyTerms(new TorsionValEnergy())[1];
		TetherEnergy tetherTerm = (TetherEnergy) energy.getEnergyTerm(new TetherEnergy());


		// Looping on all the models
		String[] models = File2StringArray.f2a(loopsFileName);
		String[] dataFile = File2StringArray.f2a(dataFileName);
		for (int i=0 ; i<models.length ; i++) {
			try {	
				// Reading the model
				AtomList loop = (new AtomList(models[i])).backbone();
				for (int c=0; c<loop.size() ; c++) {
					Atom atom = loop.atomAt(c);
					model.atoms().findAtomInList(atom.name(), atom.residueNumber()).setXYZ(atom.x(), atom.y(), atom.z());
				}
				
				// Setting the minimization energy
				tetherTerm.updatePegsToCurrentPosition();
				String[] phipsiData = File2StringArray.f2a(dataFile[i]);
				for (int c=0 ; c<phipsiData.length ; c++) {
					StringTokenizer st = new StringTokenizer(phipsiData[c]);
					int residueNumber = Integer.parseInt(st.nextToken());
					double phi = Double.parseDouble(st.nextToken());
					double psi = Double.parseDouble(st.nextToken());
					if ((residueNumber!=dontTakeTorsion1) && (residueNumber!=dontTakeTorsion2) &&
							(torsionTerm1.getTorsionEnergyElement(residueNumber, "PHI")!=null) &&
							(torsionTerm1.getTorsionEnergyElement(residueNumber, "PSI")!=null)) {  // The null checking is for minimizing terminal loops 
						torsionTerm1.getTorsionEnergyElement(residueNumber, "PHI").setTarget(phi);
						torsionTerm1.getTorsionEnergyElement(residueNumber, "PSI").setTarget(psi);
						torsionTerm2.getTorsionEnergyElement(residueNumber, "PHI").setTarget(phi);
						torsionTerm2.getTorsionEnergyElement(residueNumber, "PSI").setTarget(psi);
					}
					else {
						System.out.println("Not constraining the {phi,psi} of residue: "+residueNumber);
					}
				}

				// Energy and RMS - Before minimization
				System.out.println("999999 " + i + " 0000 " + models[i]);
				System.out.println("999999 " + i + " 1111 " + calcRMS(model, reference, resStart, resEnd) + " " + (-1));	 
				System.out.println("999999 " + i + " 2222 " + " 0 0 0");

				// Minimization
				torsionTerm1.off();
				torsionTerm2.off();
				hbTerm.off();
				minimizer = new SteepestDecent(energy, 0.1, 100, 50);
				System.out.println(minimizer.minimize());
				hbTerm.on();
				torsionTerm2.on();
				minimizer = new LBFGS(energy, 0.1, 1000, 200);
				System.out.println(minimizer.minimize());
				torsionTerm1.on();
				torsionTerm2.off();
				minimizer = new LBFGS(energy, 0.1, 3000, 200);
				System.out.println(minimizer.minimize());

				// Writing the minimized loop
				BufferedWriter bw = new BufferedWriter(new FileWriter(models[i]+ ".min"));
				for (int c=resStart; c<=resEnd ; c++)
					for (int cc=0; cc<model.residue(c).atoms().size() ; cc++)
						bw.write(model.residue(c).atoms().atomAt(cc) + "\n");
				bw.close();
				
				/* Final energies of the model */ 
				System.out.println("999999 " + i + " 3333 " + calcRMS(model, reference, resStart, resEnd) + " " +
						(-1));	 
				System.out.println("999999 " + i + " 4444 0 0 0");

			}
			catch (Exception e) {
				System.out.println("Minimization was not successful.\n");
				e.printStackTrace();
				System.out.println("\nContinueing to next model.\n");
			}
		}

	} // Of main



	private static double calcRMS(Protein prot, Protein ref, int start, int end) {
		double totRms = 0.0;
		for (int c=start; c<=end ; c++) {
			Atom atom = prot.residue(c).atoms().findAtomInList("N",c);
			Atom atomr = ref.residue(c).atoms().findAtomInList("N",c);
			totRms += (atom.x() - atomr.x())*
			(atom.x() - atomr.x()) +
			(atom.y() - atomr.y())*
			(atom.y() - atomr.y()) + 
			(atom.z() - atomr.z())*
			(atom.z() - atomr.z());
			atom = prot.residue(c).atoms().findAtomInList("CA",c);
			atomr = ref.residue(c).atoms().findAtomInList("CA",c);
			totRms += (atom.x() - atomr.x())*
			(atom.x() - atomr.x()) +
			(atom.y() - atomr.y())*
			(atom.y() - atomr.y()) + 
			(atom.z() - atomr.z())*
			(atom.z() - atomr.z());
			atom = prot.residue(c).atoms().findAtomInList("C",c);
			atomr = ref.residue(c).atoms().findAtomInList("C",c);
			totRms += (atom.x() - atomr.x())*
			(atom.x() - atomr.x()) +
			(atom.y() - atomr.y())*
			(atom.y() - atomr.y()) + 
			(atom.z() - atomr.z())*
			(atom.z() - atomr.z());
			atom = prot.residue(c).atoms().findAtomInList("O",c);
			atomr = ref.residue(c).atoms().findAtomInList("O",c);
			totRms += (atom.x() - atomr.x())*
			(atom.x() - atomr.x()) +
			(atom.y() - atomr.y())*
			(atom.y() - atomr.y()) + 
			(atom.z() - atomr.z())*
			(atom.z() - atomr.z());
		}
		return Math.sqrt(totRms/(4*(end-start+1)));
	}

	private static double calcRMSallHeavyAtoms(Protein prot, Protein ref, int start, int end) {
		double totRms = 0.0;
		int ntot = 0;
		for (int c=start; c<=end ; c++) {
			for (int d=0; d<prot.residue(c).atoms().size() ; d++) 
				if (!prot.residue(c).atoms().atomAt(d).isHydrogen) {
					Atom atom = prot.residue(c).atoms().atomAt(d);
					Atom atomr = ref.residue(c).atoms().findAtomInList(prot.residue(c).atoms().atomAt(d).name(),c);
					if (atomr!=null) {
						totRms += (atom.x() - atomr.x())*(atom.x() - atomr.x()) + 
						(atom.y() - atomr.y())*(atom.y() - atomr.y()) +
						(atom.z() - atomr.z())*(atom.z() - atomr.z());
						ntot++;
					}
				}
		}
		return Math.sqrt(totRms/ntot);
	}





	/** ================================= init =========================================
	 *
	 *A static function for parsing of the command line arguments and assigning the 
	 *variables commandsFileName, modelFileName and randomNumberSeed with the right inputs. Note that this
	 *static method is using parsing functions such as getOrderedArguments that are defined in MeshiProgram
	 *that MinimizeProtein inherits.
	 **/

	private static void init(String[] args) {

		/**** NOTE *** the next two lines. Because of a BUG in the Java VM, the 
		 * interfaces "Residues" and "AtomTypes" are not loaded automatically when MinimizeProtein initialize. 
		 * For this purpose these two lines are crucial wherever these two interfaces are implemented. The user might 
		 * rightfully feel that these two lines are "black magic" programming, but happily to our knowledge this is 
		 * the only bizarre phenomenon we are aware of in meshi.
		 **/
		int zvl = ALA; // force the reading of "meshi.parameters.Residues"
		zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"


		String line;
		String errorMessage = ("\n                  ******************\n"+
				"Usage java -Xmx300m MinimizeBatchOfCoarseLoops <commands file name> <file with list of pdbs> <file with list of phipsi data> <ref pdb file name> " +
				"<loop starting resisue> <loop ending residue> <Wev> <Whb> <Wtorval> <Wtether> <output PDB extension string>\n"+
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));
		commandsFileName = getOrderedArgument(args);
		if (commandsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# commandsFileName = "+commandsFileName);

		commands = new CommandList(commandsFileName);

		refFileName = getOrderedArgument(args);
		if (refFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# reference file name is "+refFileName);

		modelFileName = getOrderedArgument(args);
		if (modelFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# initial model file name is "+modelFileName);

		loopsFileName = getOrderedArgument(args);
		if (loopsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# initial loops file name is "+loopsFileName);

		dataFileName = getOrderedArgument(args);
		if (dataFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# data file name is "+dataFileName);

		noTorsionFileName = getOrderedArgument(args);
		if (noTorsionFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# The file with the closure torsions: "+noTorsionFileName);

		initRandom(999);

		String tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		resStart = (new Integer(tmpString)).intValue();
		System.out.println("# Starting residue is " + resStart);

		tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		resEnd = (new Integer(tmpString)).intValue();
		System.out.println("# Ending residue is " + resEnd);
	}
}
