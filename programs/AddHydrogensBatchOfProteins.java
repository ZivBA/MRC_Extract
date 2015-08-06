package programs;

import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
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

/**
 *<pre>
 * This program will minimize a batch of proteins given in a list, and will compare the results to a reference.
 * 
 * Unix usage:
 *     java -Xmx300m AddHydrogensBatchOfProteins <commands file name> <file with list of pdbs> 
 *
 * <commands file name> - A text file containing the different flags and parameters required for 
 *                        the run.
 * <file with list of pdbs> - These will be minimized
 *
 *
 **/

public class AddHydrogensBatchOfProteins extends MeshiProgram implements Residues, AtomTypes{ 

    private static CommandList commands; 
    private static String commandsFileName = null;
    private static String modelsFileName = null;  

 
    public static void main(String[] args) throws MinimizerException, LineSearchException{
	init(args); 
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
		};	
	
	// Looping on all the models
	String[] models = File2StringArray.f2a(modelsFileName);
	for (int i=0 ; i<models.length ; i++) {
		// Reading the model
		model = new Protein(new AtomList(models[i]) , new ResidueExtendedAtoms(ADD_HYDROGENS_AND_FREEZE));
		for (int cc=0 ; cc<model.atoms().size() ; cc++)
				model.atoms().atomAt(cc).setChain("A");
		
		// Energy and RMS - Before minimization
	  	model.defrost();
	  	model.atoms().filter(new AtomList.NonHydrogen()).freeze();
		distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
		energy = new TotalEnergy(model, distanceMatrix, energyCreators, commands);
		minimizer = new LBFGS(energy, 0.05, 500, 100);
		System.out.println(minimizer.minimize());
		try {
			model.atoms().print(new MeshiWriter(models[i]+".withH.pdb"));
		}
		catch (Exception e) {
			System.out.print("\nThere was a problem writing the output:\n" + e + "\n\nContinueing...\n\n");
		}
	}
	
	} // Of main

    
    
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
        }
}
