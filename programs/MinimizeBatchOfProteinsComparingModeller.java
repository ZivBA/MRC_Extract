package programs;

import meshi.applications.prediction.GDTcalculator;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.compositeTorsions.compositePropensity2DwithPP.CompositePropensity2DCreator;
import meshi.energy.compositeTorsions.ramachandranAndChi1.RamachandranAndChi1PartialResiduesCreator;
import meshi.energy.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.linearRG.LinearRgCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.softExcludedVol.SoftExcludedVolCreator;
import meshi.energy.solvate.SolvateCreatorHBforMinimization;
import meshi.energy.solvate.SolvateEnergy;
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

public class MinimizeBatchOfProteinsComparingModeller extends MeshiProgram implements Residues, AtomTypes{ 

	private static CommandList commands; 
	private static String commandsFileName = null;
	private static String modelsFileName = null;  
	private static String refFileName = null;  
	private static double Wrg = 0.0;  
	private static double Wev = 3.0;  
	private static double Wsolv = 0.5;  
	private static double Whb = 1.0;  
	private static double Wprop = 1.0;  
	private static double Wramach = 0.1;  


	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		Protein reference = null;
		Protein modelFull = null;
		Protein model = null;
		DistanceMatrix distanceMatrix = null;
		TotalEnergy energy = null;
		SolvateEnergy solvTerm = null;
		Minimizer minimizer = null;

		// The creators for the terms, and the rotamer lib
		EnergyCreator[] energyCreators = {  
				new BondCreator(),
				new AngleCreator(),
				new PlaneCreator(),
				new OutOfPlaneCreator(),
				new LinearRgCreator(Wrg),
				new SoftExcludedVolCreator(Wev , 3 , 1.0),
				new SolvateCreatorHBforMinimization(Wsolv,Wsolv,Wsolv,Whb),
				new CompositePropensity2DCreator(Wprop),
				new RamachandranSidechainEnergyCreator(Wramach),
				new RamachandranAndChi1PartialResiduesCreator(Wramach),
				new TetherCreator(0.5,new AtomList.ClassNCaCFilter())				
		};	

		EnergyCreator[] energyCreatorsForSCMODminimzation = {  
				new BondCreator(),
				new AngleCreator(),
				new PlaneCreator(),
				new OutOfPlaneCreator(),
				new SoftExcludedVolCreator(1.0 , 11 , 1.0),
				new SolvateCreatorHBforMinimization(0.01,0.01,0.01,1.0),
				new CompositePropensity2DCreator(0.01),
				new RamachandranSidechainEnergyCreator(1.0),
				new RamachandranAndChi1PartialResiduesCreator(1.0),
				new TetherCreator(0.5,new AtomList.ClassNCaCFilter())
		};	
		DunbrackLib lib = new DunbrackLib(commands,1.0,30);

		// Loading the reference	
		reference = new Protein(refFileName, new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));

		// Looping on all the models
		String[] models = File2StringArray.f2a(modelsFileName);
		for (int i=0 ; i<models.length ; i++) { 		
			// Reading the model
			System.out.println("Trying to minimize: " + models[i]);
			model = buildModel(reference,models[i]);
			if (model!=null) { 		
				try {
					// Minimizing and finding the phi,psi for SCMOD
					modelFull = new Protein(models[i], new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
					distanceMatrix = new DistanceMatrix(modelFull.atoms(), 5.5, 2.0, 4);  
					energy = new TotalEnergy(modelFull, distanceMatrix, energyCreatorsForSCMODminimzation, commands);
					minimizer = new LBFGS(energy, 0.05, 500, 100);
					System.out.println(minimizer.minimize());
					double[][] pp = RotamericTools.phipsi(modelFull, distanceMatrix);
					modelFull = new Protein(models[i], new ResidueExtendedAtoms(ADD_ATOMS));
					SCMOD.scmod(commands, lib, modelFull, 3, 0.0, pp);
					
					
					// Energy and RMS - Before minimization
					// ------------------------------------
					model = new Protein(getMatchingAtoms(reference, modelFull.atoms()),
							            new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
					System.out.println("999999 " + i + " 0000 " + models[i]);
					System.out.println("999999 " + i + " 1111 " + reference.atoms().CAFilter().getRms(model.atoms().CAFilter()) + 
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),0.5,1.0,2.0,4.0) +	 
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),1.0,2.0,4.0,8.0) + " " +	
							0.0); //reference.atoms().noOXTFilter().filter(new AtomList.NonHydrogen()).getRms(model.atoms().noOXTFilter().filter(new AtomList.NonHydrogen())));				
					/* Initial energies of the model */ 
					distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
					energy = new TotalEnergy(model, distanceMatrix, energyCreators, commands);
					energy.evaluate();
					solvTerm = (SolvateEnergy) energy.getEnergyTerm(new SolvateEnergy());
					System.out.println("999999 " + i + " 2222 " + energy.report(2) + " " + 
							solvTerm.evaluate(false,1.0,0.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,1.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,1.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,0.0,1.0));


					// Energy and RMS - After minimization
					// ------------------------------------
					distanceMatrix = new DistanceMatrix(modelFull.atoms(), 5.5, 2.0, 4);  
					energy = new TotalEnergy(modelFull, distanceMatrix, energyCreators, commands);
					minimizer = new LBFGS(energy, 0.05, 20000, 100);
					System.out.println(minimizer.minimize());

					model = new Protein(getMatchingAtoms(reference, modelFull.atoms()),
				            new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));					
					System.out.println("999999 " + i + " 3333 " + reference.atoms().CAFilter().getRms(model.atoms().CAFilter()) + 
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),0.5,1.0,2.0,4.0) +	 
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),1.0,2.0,4.0,8.0) + " " +	
							0.0); //reference.atoms().noOXTFilter().filter(new AtomList.NonHydrogen()).getRms(model.atoms().noOXTFilter().filter(new AtomList.NonHydrogen())));				
					/* Final energies of the model */ 
					distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
					energy = new TotalEnergy(model, distanceMatrix, energyCreators, commands);
					energy.evaluate();
					solvTerm = (SolvateEnergy) energy.getEnergyTerm(new SolvateEnergy());
					System.out.println("999999 " + i + " 4444 " + energy.report(2) + " " + 
							solvTerm.evaluate(false,1.0,0.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,1.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,1.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,0.0,1.0));
					try {
						model.atoms().print(new MeshiWriter(models[i]+".min.pdb"));
					}
					catch (Exception e) {
						System.out.print("\nThere was a problem writing the output:\n" + e + "\n\nContinueing...\n\n");
					}
				}
				catch (Exception e) {
					System.out.println("SKIPPING: minimization was not successful: " + models[i]);
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

	private static Protein buildModel(Protein reference, String modelName) {
		Protein tmpProt=null;
		AtomList tmpList = getMatchingAtoms(reference,	new AtomList(modelName));
		if (tmpList!=null) {
			tmpProt = new Protein(tmpList, new ResidueExtendedAtoms(ADD_ATOMS));
			tmpProt.defrost();
			for (int cc=0 ; cc<tmpProt.atoms().size() ; cc++)
				tmpProt.atoms().atomAt(cc).setChain("A");
			return tmpProt;
		}
		return tmpProt;
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
		Wrg = (new Double(tmpString)).doubleValue();
		System.out.println("# Wrg is " + Wrg);

		tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		Wev = (new Double(tmpString)).doubleValue();
		System.out.println("# Wev is " + Wev);

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
	}
}
