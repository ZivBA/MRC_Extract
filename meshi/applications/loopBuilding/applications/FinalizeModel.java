package meshi.applications.loopBuilding.applications;

import meshi.applications.corpus.Corpus;
import meshi.applications.prediction.GDTcalculator;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.LennardJones.LennardJonesCreator;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_Minimization_Creator;
import meshi.energy.tether.TetherCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.rotamers.DunbrackLib;
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
import meshi.util.file.MeshiWriter;
import programs.PutHydrogens;
import programs.SCMOD;

public class FinalizeModel extends MeshiProgram implements Residues, AtomTypes {

    private static CommandList commands; 

	private static String fileNameCorpus = null;

	private static String referenceFileName = null;
		
	private static String modelFileName = null;

	private static String outFileName = null;
	
	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		
		Protein query = null;
		Protein reference = null;
		DistanceMatrix dm = null;
		Corpus corpus = new Corpus(fileNameCorpus);
		DunbrackLib lib  = new DunbrackLib(commands, 1.0 , 90);
		
		// Loading protein and reference
		query = new Protein(new AtomList(modelFileName), new ResidueExtendedAtoms(ADD_ATOMS));
		for (int cc=0 ; cc<query.atoms().size() ; cc++)
			query.atoms().atomAt(cc).setChain(" ");
		reference = new Protein(new AtomList(referenceFileName), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		System.out.println("999999 0 1111 " + reference.atoms().CAFilter().getRms(getMatchingAtoms(reference,query.atoms()).CAFilter()) + 
				" " + GDTcalculator.gdt(reference.atoms(),query.atoms(),0.5,1.0,2.0,4.0) +	
				" " + GDTcalculator.gdt(reference.atoms(),query.atoms(),1.0,2.0,4.0,8.0) + " -1 " +
				getPercentageOfMatchingAtoms(reference.atoms(), query.atoms()));
		
		// Correcting rotten residues, and SCMODing
		System.out.println("Fixing rotten residues...");
		PutHydrogens.adjustHydrogens(commands, query);
		query.defrost();
		System.out.println("Sidechain modeling...");
		SCMOD.scmod(commands, lib, query, 2);

		// Minimizing 
		EnergyCreator[] energyCreators = {  
				new BondCreator(1.0),
				new AngleCreator(1.0),
				new PlaneCreator(10.0),
				new OutOfPlaneCreator(1.0),
				new LennardJonesCreator(0.4),
				new SimpleHydrogenBond_Dahiyat_Minimization_Creator(1.0),
				new RamachandranSidechainEnergyCreator(1.0),
				new TetherCreator(5.0, new AtomList.ClassCaCbFilter()) 
		};	
		dm = new DistanceMatrix(query.atoms(), 5.5, 2.0, 4);
//		((TetherCreator) energyCreators[7]).takePegFrom("constraints.pdb");		
		TotalEnergy energy = new TotalEnergy(query, dm, energyCreators, commands);
		Minimizer minimizer = new LBFGS(energy, 0.05, 10000 , 100);
		System.out.println(minimizer.minimize());
		System.out.println("999999 0 3333 " + reference.atoms().CAFilter().getRms(getMatchingAtoms(reference,query.atoms()).CAFilter()) + 
				" " + GDTcalculator.gdt(reference.atoms(),query.atoms(),0.5,1.0,2.0,4.0) +	
				" " + GDTcalculator.gdt(reference.atoms(),query.atoms(),1.0,2.0,4.0,8.0) + " -1 " +
				getPercentageOfMatchingAtoms(reference.atoms(), query.atoms()));
		
		// Writing to file
		try {
			query.atoms().print(new MeshiWriter(outFileName));
		}
		catch (Exception e) {
			System.out.print("\nThere was a problem writing the output:\n" + e + "\n\nContinueing...\n\n");
		}		
	}

	
	private static AtomList getMatchingAtoms(Protein refProt , AtomList protAtoms) {
		AtomList result = new AtomList();
		for (int c=0 ; c<refProt.residues().size() ; c++) 
			if (refProt.residues().residueAt(c).ca()!=null)
				if (protAtoms.findAtomInList("CA" , refProt.residues().residueAt(c).ca().residueNumber())==null) {
					System.out.println("MISMATCH: this atom is not found in model: " 
							+ refProt.residues().residueAt(c).ca());
					return null;
				}
		for (int c=0 ; c<protAtoms.size() ; c++)
			if (refProt.atoms().findAtomInList("CA" , protAtoms.atomAt(c).residueNumber())!=null)
				result.add(protAtoms.atomAt(c));

		return result;
	}

	private static double getPercentageOfMatchingAtoms(AtomList refList , AtomList modelList) {
		int n = 0;
		
		for (int c=0 ; c<modelList.size() ; c++) {
			if (modelList.atomAt(c).name().equals("CA"))
				if (refList.findAtomInList("CA", modelList.atomAt(c).residueNumber()) != null)
					n++;
		}

		return n*100.0/(refList.CAFilter().size());
	}
	
	

	/** ================================= init =========================================
	 *
	 *A static function for parsing of the command line arguments and assigning the 
	 *variables commandsFileName, modelFileName and randomNumberSeed with the right inputs. Note that this
	 *static method is using parsing functions such as getOrderedArguments that are defined in MeshiProgram
	 *that MinimizeProtein inherits.
	 **/

	protected static void init(String[] args) {

		/**** NOTE *** the next two lines. Because of a BUG in the Java VM, the 
		 * interfaces "Residues" and "AtomTypes" are not loaded automatically when MinimizeProtein initialize. 
		 * For this purpose these two lines are crucial wherever these two interfaces are implemented. The user might 
		 * rightfully feel that these two lines are "black magic" programming, but happily to our knowledge this is 
		 * the only bizarre phenomenon we are aware of in meshi.
		 **/
		int zvl = ALA; // force the reading of "meshi.parameters.Residues"
		zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"


		String errorMessage = ("\n                  ******************\n"+
				"Usage java -Xmx600m FinalizeCASP8Models <commands file name> <corpus filename> <alignment file name> <model file> <output file>\n"+
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));

		String commandsFileName = getOrderedArgument(args);
		if (commandsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# commandsFileName = "+commandsFileName);

		commands = new CommandList(commandsFileName);

		fileNameCorpus = getOrderedArgument(args);
		if (fileNameCorpus == null) throw new RuntimeException(errorMessage);
		System.out.println("# Corpus file name is "+fileNameCorpus);

		referenceFileName = getOrderedArgument(args);
		if (referenceFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Reference: "+referenceFileName);

		modelFileName = getOrderedArgument(args);
		if (modelFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Completeing: "+modelFileName);
		
		outFileName = getOrderedArgument(args);
		if (outFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Output to: "+outFileName);

		initRandom(444);
	}	

} // Of FinalizeCASP8Models
