package programs;

import meshi.applications.prediction.GDTcalculator;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.excludedVolumeImprovedDistance.EVenergyCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.solvate.SolvateCreatorRegularHB;
import meshi.energy.solvate.SolvateEnergy;
import meshi.geometry.DistanceMatrix;
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
import meshi.util.file.MeshiWriter;
import meshi.util.rotamericTools.RotamericTools;

/**
 *<pre>
 * This program will minimize a batch of proteins given in a list, and will compare the results to a reference.
 * 
 * Unix usage:
 *     java -Xmx300m MinimizeBatchOfProteins <commands file name> <file with list of pdbs> <ref pdb file name> <Wrg> <Wev> <Wsolv> <Whb> <Wprop> <Wramach>
 *
 * <commands file name> - A text file containing the different flags and parameters required for 
 *                        the run.
 * <file with list of pdbs> - These will be minimized
 *
 * <ref pdb file name> - GDT and RMS will be calculated to this structure.
 * 
 * <W...> - The weights
 *
 **/

public class PerturbNative extends MeshiProgram implements Residues, AtomTypes{ 

	private static CommandList commands; 
	private static String commandsFileName = null;
	private static String modelsFileName = null;  
	private static String refFileName = null;  
	private static double Wpol = 0.0;  
	private static double Wev = 3.0;  
	private static double Wsolv = 0.5;  
	private static double Whb = 1.0;  
	private static double Wprop = 1.0;  
	private static double Wramach = 0.1;  
	private static double frac = 1.0;  
	private static double Wtether = 1.0;  


	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		Protein reference = null;
		Protein model = null;
		DistanceMatrix distanceMatrix = null;
		TotalEnergy energy = null;
		SolvateEnergy solvTerm = null;
		Minimizer minimizer = null;		

		//Corpus corpus = new Corpus("finalCorpus.txt");
		//DunbrackLib lib  =  new DunbrackLib(commands, 1.0 , 30);


/*		// The creators for the terms
		EnergyCreator[] energyCreators1 = {  
				new BondCreator(),
				new AngleCreator(),
				new PlaneCreator(),
				new OutOfPlaneCreator(),
				new SoftExcludedVolCreator(Wev , 2 , frac),   
				new SolvateCreatorRegularHBjustBB(Wpol,Wsolv,Wpol,Whb),
				new RamachandranSidechainEnergyCreator(Wramach),
				new TetherCreator(5.0, new AtomList.ClassCaCbFilter())
		};	
*/

		// The creators for the terms
		EnergyCreator[] energyCreators2 = { 
				new BondCreator(),
				new AngleCreator(),
				new PlaneCreator(1.0),
				new EVenergyCreator(Wev,0,1.0),
//				new LennardJonesCreator(Wev , 3 , frac),
				new SolvateCreatorRegularHB(0.00001,Wsolv,0.00001,Whb),
//				new SolvateCreatorRegularHB(0.00001,Wsolv,0.00001,0.00001),
//				new SolvateCreatorRegularHBjustBB(0.00001,0.00001,0.00001,Whb),
//				new SolvateCreatorRegularHBnoBB(0.00001,0.00001,0.00001,Whb*Wtether),
				new RamachandranSidechainEnergyCreator(Wramach),
				//new TetherCreator(0.0, new KolDichfin())
				//new LinearRgCreator(10.0)
		};	

		// Loading the reference	
		reference = new Protein(refFileName, new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		RotamericTools.correctNomenclature(reference);

		// Looping on all the models
		String[] models = File2StringArray.f2a(modelsFileName);
		for (int i=0 ; i<models.length ; i++) { 		
			// Reading the model
			System.out.println("Trying to minimize: " + models[i]);
			//AtomList tmpList = getMatchingAtoms(reference,	new AtomList(models[i]));
			//if (tmpList!=null) {
			if (true) {
				//model = new Protein(tmpList, new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
				model = new Protein(new AtomList(models[i]), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
				for (int cc=0 ; cc<model.atoms().size() ; cc++)
					model.atoms().atomAt(cc).setChain("A");
				model.defrost();

				// Energy and RMS - Before minimization
				System.out.println("999999 " + i + " 0000 " + models[i]);
				System.out.println("999999 " + i + " 1111 " + reference.atoms().CAFilter().getRms(getMatchingAtoms(reference,model.atoms()).CAFilter()) + 
						" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),0.5,1.0,2.0,4.0) +	
						" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),1.0,2.0,4.0,8.0) + " " +	
						0.0); //reference.atoms().noOXTFilter().filter(new AtomList.NonHydrogen()).getRms(model.atoms().noOXTFilter().filter(new AtomList.NonHydrogen())));				
				distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
				energy = new TotalEnergy(model, distanceMatrix, energyCreators2, commands);
				//energy = new TotalEnergyTorsionSpace(model, distanceMatrix, energyCreators2, commands);
				energy.evaluate();
				solvTerm = (SolvateEnergy) energy.getEnergyTerm(new SolvateEnergy());
				System.out.println("999999 " + i + " 2222 " + energy.report(2) + " " + 
						solvTerm.evaluate(false,1.0,0.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,1.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,1.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,0.0,1.0));					

				minimizer = new LBFGS(energy, 0.1, 100000, 100);
				try {
				System.out.println(minimizer.minimize());
				}
				catch (Exception e) {
					System.out.print("\nThere was a problem in minimizing.\n" + e + "\n\nContinueing...\n\n");
				}
/*				System.out.println("999999 " + i + " 5555 " + reference.atoms().CAFilter().getRms(getMatchingAtoms(reference,model.atoms()).CAFilter()) + 
						" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),0.5,1.0,2.0,4.0) +	 
						" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),1.0,2.0,4.0,8.0) + " " +	
						0.0); //reference.atoms().noOXTFilter().filter(new AtomList.NonHydrogen()).getRms(model.atoms().noOXTFilter().filter(new AtomList.NonHydrogen())));				

				System.out.println("999999 " + i + " 6666 " + energy.report(2) + " " + 
						solvTerm.evaluate(false,1.0,0.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,1.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,1.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,0.0,1.0));
				
				//energy.getEnergyTerm(new TetherEnergy()).off();				
				minimizer = new LBFGS(energy, 0.001, 30000, 100);
				System.out.println(minimizer.minimize());
*/				
				System.out.println("999999 " + i + " 3333 " + reference.atoms().CAFilter().getRms(getMatchingAtoms(reference,model.atoms()).CAFilter()) + 
						" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),0.5,1.0,2.0,4.0) +	 
						" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),1.0,2.0,4.0,8.0) + " " +	
						0.0); //reference.atoms().noOXTFilter().filter(new AtomList.NonHydrogen()).getRms(model.atoms().noOXTFilter().filter(new AtomList.NonHydrogen())));				
				System.out.println("999999 " + i + " 4444 " + energy.report(2) + " " + 
						solvTerm.evaluate(false,1.0,0.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,1.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,1.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,0.0,1.0));
				try {
					model.atoms().print(new MeshiWriter(models[i]+".min.pdb"));
				}
				catch (Exception e) {
					System.out.print("\nThere was a problem writing the output:\n" + e + "\n\nContinueing...\n\n");
				}
			}
		}
	} // Of main




private static AtomList getMatchingAtoms(Protein refProt , AtomList protAtoms) {
	AtomList result = new AtomList();
	for (int c=0 ; c<refProt.residues().size() ; c++) 
		if (refProt.residues().residueAt(c).ca()!=null)
			if (protAtoms.findAtomInList("CA" , refProt.residues().residueAt(c).ca().residueNumber())==null) {
				System.out.println("SKIPPING: this atom not found in model: " 
						+ refProt.residues().residueAt(c).ca());
				return null;
			}
	for (int c=0 ; c<protAtoms.size() ; c++)
		if (refProt.atoms().findAtomInList("CA" , protAtoms.atomAt(c).residueNumber())!=null)
			result.add(protAtoms.atomAt(c));

	return result;
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
			"Usage java -Xmx300m MinimizeBatchOfProteins <commands file name> <file with list of pdbs> <ref pdb file name> " +
			"<Wrg>  <Wev> <Wsolv> <Whb> <Wprop> <Wramach> \n"+
	"                    ******************\n");

	if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));
	commandsFileName = getOrderedArgument(args);
	if (commandsFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# commandsFileName = "+commandsFileName);

	commands = new CommandList(commandsFileName);

	modelsFileName = getOrderedArgument(args);
	if (modelsFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# initial model file name is "+modelsFileName);

	refFileName = getOrderedArgument(args);
	if (refFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# reference file name is "+refFileName);

	initRandom(999);

	String tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	Wev = (new Double(tmpString)).doubleValue();
	System.out.println("# Wev is " + Wev);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	Wpol = (new Double(tmpString)).doubleValue();
	System.out.println("# Wpol is " + Wpol);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	Wsolv = (new Double(tmpString)).doubleValue();
	System.out.println("# Wsolv is " + Wsolv);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	Whb = (new Double(tmpString)).doubleValue();
	System.out.println("# Whb is " + Whb);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	Wprop = (new Double(tmpString)).doubleValue();
	System.out.println("# Wprop is " + Wprop);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	Wramach = (new Double(tmpString)).doubleValue();
	System.out.println("# Wramach is " + Wramach);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	Wtether = (new Double(tmpString)).doubleValue();
	System.out.println("# Wtether is " + Wtether);

	tmpString = getOrderedArgument(args);
	if (tmpString!= null) {
		frac = (new Double(tmpString)).doubleValue();
		System.out.println("# EV frac " + frac);
	}
}
}
