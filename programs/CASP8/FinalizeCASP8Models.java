package programs.CASP8;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.excludedVolumeImprovedDistance.EVenergyCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.solvateNew.SolvateCreatorRegularHB;
import meshi.energy.tether.TetherCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
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
import meshi.util.file.MeshiWriter;
import meshi.util.rotamericTools.RotamericTools;
import programs.PutHydrogens;
import programs.SCMOD;

public class FinalizeCASP8Models extends MeshiProgram implements Residues, AtomTypes {

    private static CommandList commands; 

	private static String alignmentFileName = null;  

	private static String fileNameCorpus = null;
	
	private static String modelFileName = null;

	private static String outFileName = null;
	
	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		
		Protein query = null;
		String queryAlignment = "";
		int resNumCounterQ=0;
		DistanceMatrix dm = null;

		
		// Reading the alignment file
		String[] alignmentStrings = File2StringArray.f2a(alignmentFileName);
		queryAlignment = alignmentStrings[2].trim().replace("-", "");
		System.out.println("SEQ: " + queryAlignment);

		
		// Adding missing residues
		Protein oldModel = new ExtendedAtomsProtein(modelFileName,DO_NOT_ADD_ATOMS);
		AtomList largerList = new AtomList();
		resNumCounterQ = 0;
		Atom lastAdded=oldModel.residue(oldModel.firstResidue()).ca();
		for (int position=0; position<queryAlignment.length() ; position++) {
			resNumCounterQ++;
			if ((oldModel.residue(resNumCounterQ)!=null) &&
					!oldModel.residue(resNumCounterQ).dummy()) { // Coping atoms from existing model
				for (int atomCounter=0 ; atomCounter<oldModel.residue(resNumCounterQ).atoms().size() ; atomCounter++) { 
					largerList.add(new Atom(oldModel.residue(resNumCounterQ).atoms().atomAt(atomCounter)));		
					lastAdded = oldModel.residue(resNumCounterQ).atoms().atomAt(atomCounter);
				}
			}
			else {	
				System.out.println("Added residue: "+Residue.one2three(queryAlignment.charAt(position))+" " +
						resNumCounterQ);
				largerList.add(new Atom(lastAdded.x(),lastAdded.y(),lastAdded.z(),1.0,"CA",
						Residue.one2three(queryAlignment.charAt(position)),
						resNumCounterQ,-1));
			}
		}
		
		
		query = new Protein(largerList, new ResidueExtendedAtoms(ADD_ATOMS));
		for (int cc=0 ; cc<query.atoms().size() ; cc++)
			query.atoms().atomAt(cc).setChain(" ");
		PutHydrogens.adjustHydrogens(commands, query);		 

//		EnergyCreator[] energyCreators1 = {  
//				new BondCreator(),
//				new AngleCreator(),
//				new PlaneCreator(),
//				new OutOfPlaneCreator(),
//				new EVenergyCreator(0.05,0,1.0),
//				new TetherCreator(100.0, new AtomList.NonHydrogen()) 
//		};			
//		// Some quick run
//		query.defrost();
//		dm = new DistanceMatrix(query.atoms(), 5.5, 2.0, 4);
//		((TetherCreator) energyCreators1[5]).takePegFrom(modelFileName);		
//		TotalEnergy energy1 = new TotalEnergy(query, dm, energyCreators1, commands);
//		energy1.evaluate();
//		Minimizer minimizer1 = new LBFGS(energy1, 0.05, 251, 100);
//		System.out.println(minimizer1.minimize());
//		
//
//		DunbrackLib lib  =  new DunbrackLib(commands, 1.0 , 90);
//		SCMOD.scmod(commands, lib, query, 2, 0.0);
		// Minimizing with a strong tether
		EnergyCreator[] energyCreators = {  
				new BondCreator(),
				new AngleCreator(),
				new PlaneCreator(),
				new OutOfPlaneCreator(),
				new EVenergyCreator(0.05,0,1.0),
				new SolvateCreatorRegularHB(0.0,4.0,0.0,1.2),
				new RamachandranSidechainEnergyCreator(0.4),
				new TetherCreator(3.0, new AtomList.ClassCaCbFilter()) 
		};			
		// Some quick run
		query.defrost();		
//		for (int t=1; t<49; t++)
//			query.residue(t).freeze();
		dm = new DistanceMatrix(query.atoms(), 5.5, 2.0, 4);
		((TetherCreator) energyCreators[7]).takePegFrom(modelFileName);		
//		((TetherCreator) energyCreators[7]).takePegFrom("constraints.pdb");		
		TotalEnergy energy = new TotalEnergy(query, dm, energyCreators, commands);
		energy.evaluate();
		Minimizer minimizer = new LBFGS(energy, 0.05, 1000 , 100);
		System.out.println(minimizer.minimize());
		
		// Writing to file
		try {
			query.atoms().print(new MeshiWriter(outFileName));
		}
		catch (Exception e) {
			System.out.print("\nThere was a problem writing the output:\n" + e + "\n\nContinueing...\n\n");
		}		

		// Fixing rotten residues, and giving them SC modeling
		System.out.println("Fixing rotten residues...");
//		Corpus corpus = new Corpus(fileNameCorpus);
		DunbrackLib lib  =  new DunbrackLib(commands, 1.0 , 90);
		dm = new DistanceMatrix(query.atoms(), 5.5, 2.0, 4);
		double[][] pp = RotamericTools.phipsi(query, dm);
//		CorrectRottenRamachandrans corr = new CorrectRottenRamachandrans(commands , corpus, lib);
//		double[][] noRot = corr.detectAndCorrectRottenResidues(query);
//		// Fixing the changed side-chains
//		query.freeze();
//		for (int c=0 ; c<noRot.length; c++) {
//			pp[(int) noRot[c][0]][0] = noRot[c][1]; // new Phi
//			pp[(int) noRot[c][0]][1] = noRot[c][2]; // new Psi
//			query.residue((int) noRot[c][0]).atoms().defrost();
//		}
		System.out.println("Sidechain modeling...");
		query.defrost();
		SCMOD.scmod(commands, lib, query, 2, 0.0, pp);
		
		query.defrost();
		dm = new DistanceMatrix(query.atoms(), 5.5, 2.0, 4);
		((TetherCreator) energyCreators[7]).takePegFrom(modelFileName);		
//		((TetherCreator) energyCreators[7]).takePegFrom("constraints.pdb");		
		energy = new TotalEnergy(query, dm, energyCreators, commands);
//		energy.getEnergyTerm(new TetherEnergy()).off();
		energy.evaluate();
		minimizer = new LBFGS(energy, 0.05, 20000, 100);
		System.out.println(minimizer.minimize());
			
		// Writing to file
		try {
			query.atoms().print(new MeshiWriter(outFileName));
		}
		catch (Exception e) {
			System.out.print("\nThere was a problem writing the output:\n" + e + "\n\nContinueing...\n\n");
		}		
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

		alignmentFileName = getOrderedArgument(args);
		if (alignmentFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Alignment file name is: "+alignmentFileName);

		modelFileName = getOrderedArgument(args);
		if (modelFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Completeing: "+modelFileName);
		
		outFileName = getOrderedArgument(args);
		if (outFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Output to: "+outFileName);

		initRandom(444);
	}	

} // Of FinalizeCASP8Models
