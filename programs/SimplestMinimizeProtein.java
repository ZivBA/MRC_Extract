package programs;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.LennardJones.LennardJonesCreator;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.simpleHPterm.SimpleHPCreator;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBondEnergy;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_Minimization_NoDuplications_Creator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.optimizers.LBFGS;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.Minimizer;
import meshi.optimizers.MinimizerException;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;


public class SimplestMinimizeProtein extends MeshiProgram implements Residues, AtomTypes{ 
    private static CommandList commands; 
   
    private static String commandsFileName = null;
 
    private static String modelFileName = null;  

    private static Protein protein;

 
    public static void main(String[] args) throws MinimizerException, LineSearchException{
	init(args); 

	protein = new ExtendedAtomsProtein(modelFileName,DO_NOT_ADD_ATOMS); 
	protein.defrost();
		
	DistanceMatrix distanceMatrix = new DistanceMatrix(protein.atoms(), 5.5,  2.0,  4);  

	EnergyCreator[] energyCreators = {  
	    new BondCreator(),
	    new AngleCreator(),
	    new PlaneCreator(),
	    new OutOfPlaneCreator(),
	    new LennardJonesCreator(0.5),
		new RamachandranSidechainEnergyCreator(0.4),
	    new SimpleHydrogenBond_Dahiyat_Minimization_NoDuplications_Creator(1.0),
	    new SimpleHPCreator(0.01,0.01,4.5,4.5,true)
	};

	TotalEnergy energy = new TotalEnergy(protein, distanceMatrix, energyCreators, commands);
	
//	AbstractEnergy hp = energy.getEnergyTerm(new SimpleHP());
//	hp.evaluateAtoms();
//	protein.atoms().print();
//	System.exit(1);
	
	Minimizer minimizer = new LBFGS(energy, 0.001 , 100000 , 100 );  

    System.out.println(minimizer.minimize());
    
    protein.atoms().print();

    ((SimpleHydrogenBondEnergy) energy.getEnergyTerm(new SimpleHydrogenBondEnergy())).printHBlist();

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
			       "Usage java -Xmx300m SimplestMinimizeProtein <commands file name> <pdb file name> \n"+
			       "                    ******************\n");
			      
	if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));
	commandsFileName = getOrderedArgument(args);
	if (commandsFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# commandsFileName = "+commandsFileName);

	commands = new CommandList(commandsFileName);
	
	modelFileName = getOrderedArgument(args);
	if (modelFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# initial model file name is "+modelFileName);

	initRandom(999);
    }
}
