package meshi.applications.loopBuilding.applications;

import java.io.IOException;

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
import meshi.geometry.DistanceMatrix;
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
import meshi.util.file.MeshiWriter;
import meshi.util.filters.KolDichfin;


public class RefinePreModel extends MeshiProgram implements Residues, AtomTypes{ 

	private static CommandList commands_init; 
	private static String modelFileName_init = null;  
	private static String tempFileName_init = null;  	
	private static String outputFileName_init = null;  

	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		refine(commands_init, modelFileName_init, tempFileName_init, outputFileName_init);


	} // Of main



	public static void refine(CommandList commands, String modelFileName, String tempFileName, String outputFileName) throws MinimizerException, LineSearchException {	
		Protein model = null;
		DistanceMatrix distanceMatrix = null;
		TotalEnergy energy = null;
		Minimizer minimizer = null;


		// Solving creator
		EnergyCreator[] energyCreatorsAll = {  
				new SoftExcludedVolCreator(10.0 , 1 , 1.0),
				new TetherCreator(100, new KolDichfin())
		};	
		
		
		// The creators for the terms
		EnergyCreator[] energyCreators = {  
				new BondCreator(),
				new AngleCreator(),
				new PlaneCreator(),
				new OutOfPlaneCreator(),
				new SoftExcludedVolCreator(10.0 , 12 , 1.0),
				new SimpleHydrogenBond_Dahiyat_Minimization_BBonly_Creator(3.0),
				new RamachandranCreator(1.0),
				new TetherCreator(1.0, new AtomList.BackboneFilter())
		};	

		// Loading the reference and the background model	
		model = new Protein((new AtomList(modelFileName)).noOXTFilter().backbone(), new ResidueExtendedAtoms(ADD_HYDROGENS_AND_FREEZE));
		
		// Short minimization to solve problems
		model.defrost();
		TetherCreator tetherCreator1 = ((TetherCreator) energyCreatorsAll[1]);
		tetherCreator1.takePegFrom(tempFileName);
		distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
		energy = new TotalEnergy(model, distanceMatrix, energyCreatorsAll, commands);
		minimizer = new SteepestDecent(energy, 0.1, 100, 50);
		System.out.println(minimizer.minimize());
//		PutHydrogens.adjustHydrogens(commands, model);
		// End of short minimization		
		
		for (int cc=0 ; cc<model.atoms().size() ; cc++)
			model.atoms().atomAt(cc).setChain("A");
		model.defrost();
		TetherCreator tetherCreator = ((TetherCreator) energyCreators[7]);
		tetherCreator.takePegFrom(tempFileName);
		distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
		energy = new TotalEnergy(model, distanceMatrix, energyCreators, commands);

		// Minimizing
		minimizer = new SteepestDecent(energy, 0.1, 100, 50);
		System.out.println(minimizer.minimize());
		minimizer = new LBFGS(energy, 0.01, 5000, 200);
		System.out.println(minimizer.minimize());

		// Outputing
		try {
			model.atoms().print(new MeshiWriter(outputFileName));
		} catch (IOException e) {
			System.out.println("\n\nWriting refined model to disk was unsuccessful\n\n");
		}
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
				"Usage java -Xmx300m RefinePreModel <commands file name> <PreModelFile with all residues> <PreModelFile with gaps> <output file name>\n"+
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));

		String commandsFileName = getOrderedArgument(args);
		if (commandsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# commandsFileName = "+commandsFileName);
		commands_init = new CommandList(commandsFileName);

		modelFileName_init = getOrderedArgument(args);
		if (modelFileName_init == null) throw new RuntimeException(errorMessage);
		System.out.println("# initial model file name is "+modelFileName_init);

		tempFileName_init = getOrderedArgument(args);
		if (tempFileName_init == null) throw new RuntimeException(errorMessage);
		System.out.println("# Template like model file name is "+tempFileName_init);

		outputFileName_init = getOrderedArgument(args);
		if (outputFileName_init == null) throw new RuntimeException(errorMessage);
		System.out.println("# Output file name is "+outputFileName_init);

		initRandom(999);
	}
}
