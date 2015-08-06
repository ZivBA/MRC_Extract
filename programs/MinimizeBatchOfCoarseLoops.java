package programs;

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
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_Minimization_BBonly_Creator;
import meshi.energy.softExcludedVol.SoftExcludedVolCreator;
import meshi.energy.tether.TetherCreator;
import meshi.energy.torsionVal.TorsionValCreator;
import meshi.energy.torsionVal.TorsionValEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LBFGS;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.Minimizer;
import meshi.optimizers.MinimizerException;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;
import meshi.util.rotamericTools.RotamericTools;

/**
 *<pre>
 * This program will minimize a batch of loops given in a list, and will compare the results to a reference.
 * Only the loop residues are compared in the RMS. There is no superposition, and it is assumed the rest of
 * the protein is aligned with the reference.
 * Unix usage:
 *     java -Xmx300m MinimizeBatchOfLoops <commands file name> <file with list of pdbs> <ref pdb file name> <loop starting resisue> <loop ending residue> <Wrg> <Wev> <Wsolv> <Whb> <Wprop> <Wramach> <output PDB extension string>
 *
 * <commands file name> - A text file containing the different flags and parameters required for 
 *                        the run.
 * <file with list of pdbs> - These will be minimized
 *
 * <ref pdb file name> - GDT and RMS will be calculated to this structure.
 * 
 * <W...> - The weights
 * 
 * <output PDB extension string> - Each minimized PDB will be saved with this string folowing its original name 
 *
 **/

public class MinimizeBatchOfCoarseLoops extends MeshiProgram implements Residues, AtomTypes{ 

	private static CommandList commands; 
	private static String commandsFileName = null;
	private static String modelsFileName = null;  
	private static String dataFileName = null;  
	private static String refFileName = null;  
	private static String endString = null;
	private static int resStart = -999;
	private static int resEnd = -999;
	private static double Wev = -1;  
	private static double Whb = -1;  
	private static double Wtorval = -1;  
	private static double Wtether = -1;  


	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		DunbrackLib lib = new DunbrackLib(commands, 1.0, 1);
		Protein reference = null;
		Protein model = null;
		DistanceMatrix distanceMatrix = null;
		TotalEnergy energy = null;
		Minimizer minimizer = null;


		// The creators for the terms
		EnergyCreator[] energyCreators = {  
				new BondCreator(),
				new AngleCreator(),
				new PlaneCreator(),
				new OutOfPlaneCreator(),
				new SoftExcludedVolCreator(Wev , 12 , 1.0),
				new SimpleHydrogenBond_Dahiyat_Minimization_BBonly_Creator(Whb),
				new TorsionValCreator(Wtorval),
				new TorsionValCreator(50*Wtorval),
				new RamachandranCreator(1.0),
				new TetherCreator(Wtether, new AtomList.ClassCaCbFilter())
		};	


		// Loading the reference and the background model	
		reference = new Protein((new AtomList(refFileName)).noOXTFilter(), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		model = new Protein((new AtomList(refFileName)).noOXTFilter(), new ResidueExtendedAtoms(ADD_ATOMS));
		for (int cc=0 ; cc<model.atoms().size() ; cc++)
			model.atoms().atomAt(cc).setChain(" ");
		RotamericTools.putIntoRot1(model, new DistanceMatrix(model.atoms(), 5.5, 2.0, 4), lib);
		model.freeze();
		for (int c=resStart ; c<=resEnd ; c++)
			model.residue(c).atoms().defrost();
		distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
		energy = new TotalEnergy(model, distanceMatrix, energyCreators, commands);
		TorsionValEnergy torsionTerm1 = (TorsionValEnergy) energy.getEnergyTerms(new TorsionValEnergy())[0];
		TorsionValEnergy torsionTerm2 = (TorsionValEnergy) energy.getEnergyTerms(new TorsionValEnergy())[1];


		// Looping on all the models
		String[] models = File2StringArray.f2a(modelsFileName);
		String[] dataFile = File2StringArray.f2a(dataFileName);
		for (int i=0 ; i<models.length ; i++) {
			try {	
				// Reading the model
				AtomList loop = new AtomList(models[i]);
				for (int c=0; c<loop.size() ; c++) {
					Atom atom = loop.atomAt(c);
					model.atoms().findAtomInList(atom.name(), atom.residueNumber()).setXYZ(atom.x(), atom.y(), atom.z());
				}
				// Energy and RMS - Before minimization
				System.out.println("999999 " + i + " 0000 " + models[i]);
				System.out.println("999999 " + i + " 1111 " + calcRMS(model, reference, resStart, resEnd) + " " +
						calcRMSallHeavyAtoms(model, reference, resStart, resEnd));	 
				// Setting the minimization energy
				String[] phipsiData = File2StringArray.f2a(dataFile[i]);
				for (int c=0 ; c<phipsiData.length ; c++) {
					StringTokenizer st = new StringTokenizer(phipsiData[c]);
					int residueNumber = (new Integer(st.nextToken())).intValue();
					int resType = model.residue(residueNumber).type;
					double phi = (new Double(st.nextToken())).doubleValue();
					double psi = (new Double(st.nextToken())).doubleValue();
					torsionTerm1.getTorsionEnergyElement(residueNumber, "PHI").setTarget(phi);
					torsionTerm1.getTorsionEnergyElement(residueNumber, "PSI").setTarget(psi);
					torsionTerm2.getTorsionEnergyElement(residueNumber, "PHI").setTarget(phi);
					torsionTerm2.getTorsionEnergyElement(residueNumber, "PSI").setTarget(psi);
					for (int cc=0 ; cc<lib.getChiMax(resType) ; cc++) {
						double[] rot = lib.getRotamer(resType, phi, psi, 0);
						torsionTerm1.getTorsionEnergyElement(residueNumber, "CHI"+(cc+1)).setTarget(rot[cc]);
						torsionTerm2.getTorsionEnergyElement(residueNumber, "CHI"+(cc+1)).setTarget(rot[cc]);
					}
				}
				System.out.println("999999 " + i + " 2222 " + " 0 0 0");

				// Energy and RMS - After minimization
				torsionTerm1.off();
				torsionTerm2.on();
				minimizer = new LBFGS(energy, 0.1, 1000, 500);
				System.out.println(minimizer.minimize());
				torsionTerm1.on();
				torsionTerm2.off();
				minimizer = new LBFGS(energy, 0.1, 2000, 500);
				System.out.println(minimizer.minimize());

				System.out.println("999999 " + i + " 3333 " + calcRMS(model, reference, resStart, resEnd) + " " +
						calcRMSallHeavyAtoms(model, reference, resStart, resEnd));	 

				/* Final energies of the model */ 
				System.out.println("999999 " + i + " 4444 0 0 0");

				// Writing the miniizzed loop
				try {
					BufferedWriter bw = new BufferedWriter(new FileWriter(models[i]+ "." + endString));
					for (int c=resStart; c<=resEnd ; c++)
						for (int cc=0; cc<model.residue(c).atoms().size() ; cc++)
							bw.write(model.residue(c).atoms().atomAt(cc) + "\n");
					bw.close();
				}
				catch (Exception e) {
					System.out.print("\nThere was a problem writing the output:\n" + e + "\n\nContinueing...\n\n");
				}
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
		}
		return Math.sqrt(totRms/(3*(end-start+1)));
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

		modelsFileName = getOrderedArgument(args);
		if (modelsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# initial model file name is "+modelsFileName);

		dataFileName = getOrderedArgument(args);
		if (dataFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# data file name is "+dataFileName);

		refFileName = getOrderedArgument(args);
		if (refFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# reference file name is "+refFileName);

		initRandom(999);

		String tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		resStart = (new Integer(tmpString)).intValue();
		System.out.println("# Starting residue is " + resStart);

		tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		resEnd = (new Integer(tmpString)).intValue();
		System.out.println("# Ending residue is " + resEnd);


		tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		Wev = (new Double(tmpString)).doubleValue();
		System.out.println("# Wev is " + Wev);

		tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		Whb = (new Double(tmpString)).doubleValue();
		System.out.println("# Whb is " + Whb);

		tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		Wtorval = (new Double(tmpString)).doubleValue();
		System.out.println("# Wtorval is " + Wtorval);

		tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		Wtether = (new Double(tmpString)).doubleValue();
		System.out.println("# Wtether is " + Wtether);

		endString = getOrderedArgument(args);
		if (endString== null) throw new RuntimeException(errorMessage);
		System.out.println("# Each file name will end with: " + endString);
	}
}
