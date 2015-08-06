package programs;

import java.text.DecimalFormat;
import java.util.StringTokenizer;

import meshi.applications.corpus.PropensityMatrix;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.ROT1solvation.CentroidSolvationCreator;
import meshi.energy.ROT1solvation.ROT1SolvationEnergy;
import meshi.energy.compositeTorsions.ramachandran.RamachandranCreator;
import meshi.energy.compositeTorsions.ramachandran.RamachandranEnergy;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBondEnergy;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_HighAccuracy_BBonly_Creator;
import meshi.energy.softExcludedVol.SoftExcludedVol;
import meshi.energy.softExcludedVol.SoftExcludedVolCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.ResidueBuilder;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LineSearchException;
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

public class EvaluateBatchOfCoarseLoopsSecondAttempt extends MeshiProgram implements Residues, AtomTypes{ 

	private static CommandList commands; 
	private static String commandsFileName = null;
	private static String modelFileName = null;  
	private static String loopsFileName = null;  
	private static String dataFileName = null;  
	private static String refFileName = null;  
	private static int resStart = -999;
	private static int resEnd = -999;


	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		Protein reference = null;
		Protein model = null;
		DistanceMatrix distanceMatrix = null;
		TotalEnergy energy = null;
		DecimalFormat fmt2 = new DecimalFormat("0.###");


		EnergyCreator[] energyCreators = {  
//				new BondCreator(),
//				new AngleCreator(),
//				new PlaneCreator(),
//				new OutOfPlaneCreator(),
				new SoftExcludedVolCreator(1.0 , 12 , 1.0),
				new SimpleHydrogenBond_Dahiyat_HighAccuracy_BBonly_Creator(1.0),
				new RamachandranCreator(1.0),
				new CentroidSolvationCreator(1.0,"8.0",false)
		};	


		// Loading the reference and the background model	
//		reference = new Protein((new AtomList(refFileName)).noOXTFilter().backbone(), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
//		model = new Protein((new AtomList(modelFileName)).noOXTFilter().backbone(), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		reference = new Protein((new AtomList(refFileName)).noOXTFilter().backbone(), new ResidueExtendedAtoms(ADD_HYDROGENS_AND_FREEZE));
		PutHydrogens.adjustHydrogens(commands, reference);
		model = new Protein((new AtomList(modelFileName)).noOXTFilter().backbone(), new ResidueExtendedAtoms(ADD_HYDROGENS_AND_FREEZE));
		PutHydrogens.adjustHydrogens(commands, model);
		reference.defrost();
		model.defrost();
		for (int res=0 ; res<reference.residues().size(); res++) {
			if ((reference.residues().residueAt(res).type<20) &&
					(reference.residues().residueAt(res).type>-1)) {
				try {
					ResidueBuilder.buildCentroid(reference.residues().residueAt(res));
				}
				catch (Exception ee) {}
			}
		}		
		for (int cc=0 ; cc<model.atoms().size() ; cc++)
			model.atoms().atomAt(cc).setChain("A");
		double[][] originalPositions = new double[model.atoms().size()][3];
		for (int cc=0 ; cc<model.atoms().size() ; cc++) {
			originalPositions[cc][0] = model.atoms().atomAt(cc).x();
			originalPositions[cc][1] = model.atoms().atomAt(cc).y();
			originalPositions[cc][2] = model.atoms().atomAt(cc).z();
		}
		
		
		
		distanceMatrix = new DistanceMatrix(model.atoms(), 8.0, 1.0, 4);  
		energy = new TotalEnergy(model, distanceMatrix, energyCreators, commands);
		SimpleHydrogenBondEnergy hbEnergy = (SimpleHydrogenBondEnergy) energy.getEnergyTerm(new SimpleHydrogenBondEnergy());
		RamachandranEnergy ramachEnergy = (RamachandranEnergy) energy.getEnergyTerm(new RamachandranEnergy());
		ROT1SolvationEnergy CentTerm = (ROT1SolvationEnergy) energy.getEnergyTerm(new ROT1SolvationEnergy());
		SoftExcludedVol evTerm = (SoftExcludedVol) energy.getEnergyTerm(new SoftExcludedVol());
		PropensityMatrix propMatrix = new PropensityMatrix(commands,"allResPropensity.dat");
				

		// Looping on all the models
		String[] models = File2StringArray.f2a(loopsFileName);
		String[] dataFile = File2StringArray.f2a(dataFileName);

		for (int i=-1 ; i<models.length ; i++) {
			try {	
				// Reading the model
				if (i>-1) {
					AtomList loop = (new AtomList(models[i])).backbone();
					for (int c=0; c<loop.size() ; c++) {
						Atom atom = loop.atomAt(c);
						model.atoms().findAtomInList(atom.name(), atom.residueNumber()).setXYZ(atom.x(), atom.y(), atom.z());
					}
				}

				// Setting the minimization energy
				double clusterSize = 0;
				if (i>-1) {
					String[] phipsiData = File2StringArray.f2a(dataFile[i]);
					StringTokenizer st = new StringTokenizer(phipsiData[0]);
					st.nextToken();
					st.nextToken();
					st.nextToken();
					clusterSize = Double.parseDouble(st.nextToken());
				}

				energy.update();
				double[][] pp = RotamericTools.phipsi(model, distanceMatrix);
				double Ehb = hbEnergy.evaluate();
				double Eramach = ramachEnergy.evaluate();
				double Eev = evTerm.evaluate();
				double Eprop = calcProp(pp,propMatrix);
				for (int res=0 ; res<model.residues().size(); res++) {
					if ((model.residues().residueAt(res).type<20) &&
							(model.residues().residueAt(res).type>-1)) {
						try {
							ResidueBuilder.buildCentroid(model.residues().residueAt(res));
						}
						catch (Exception ee) {}
					}
				}		
				energy.update();
				double Ecent = CentTerm.evaluate();
				double rmsAll = calcRMSallHeavyAtoms(model, reference, resStart, resEnd);
				double rmsBB = calcRMS(model, reference, resStart, resEnd);
				for (int cc=0 ; cc<model.atoms().size() ; cc++) {
					model.atoms().atomAt(cc).setXYZ(originalPositions[cc][0], 
							originalPositions[cc][1], 
							originalPositions[cc][2]);
				}
				
				System.out.println("999999 " + i + " " + fmt2.format(rmsBB) + " " + fmt2.format(rmsAll) + " " +
						fmt2.format(Ecent) + " " + fmt2.format(Eramach) + " " + fmt2.format(Ehb) + " " +
						fmt2.format(Eev) + " " + fmt2.format(Eprop) + " " + clusterSize/*fmt2.format(Math.log(clusterSize+1))*/);
			}
			catch (Exception e) {
				System.out.println("Evaluation was not successful.\n");
				e.printStackTrace();
				System.out.println("\nContinueing to next model.\n");
			}
		}

	} // Of main


	
	private static double calcProp(double[][] pp, PropensityMatrix propMatrix) {
		double sum = 0.0;
		for (int res=resStart-1 ; res<=(resEnd+1) ; res++) {
			sum += propMatrix.propVal((int) pp[res][2], pp[res][0], pp[res][1]);			
		}
		return sum;
	}
	

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
		System.out.println("# Loops loops file name is "+loopsFileName);

		dataFileName = getOrderedArgument(args);
		if (dataFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# data file name is "+dataFileName);

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




