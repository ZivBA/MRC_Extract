package meshi.applications.corpus;

import meshi.energy.EnergyCreator;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;

public class BuildCorpus extends MeshiProgram implements Residues, 
							     AtomTypes{ /**
									 * The implemented
									 * interfaces defines the 
									 * names of atom and residue 
									 * types. 
									 **/
    /**
     * Reads, parse and stores the contents of the commands file. 
     **/ 
    private static CommandList commands; 
   
    /**
     * A string with the name of the command file.
     **/
    private static String commandsFileName = null;
 
    /**
     * A string with the name of the pdb file to minimize.
     **/
    private static String modelFileName = null;  


    public static void main(String[] args) {
	init(args); 
	EnergyCreator[] energyCreators = {  
	   // new CompositePropensity2DCreator(),
	};
	Corpus corpus = new Corpus(modelFileName,commands,energyCreators);
	corpus.writeToDisk(modelFileName + ".corpus");

	//Corpus corpus2 = new Corpus(modelFileName + ".corpus");
	//corpus2.writeToDisk(modelFileName + ".corpus2");
	//corpus.merge(corpus2);
	//corpus.writeToDisk(modelFileName + ".corpus.merge");
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
			       "Usage java -Xmx300m BuildCorpus <commands file name> <pdb file name> seed\n"+
			       "                    ******************\n");
			      
	if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));
	commandsFileName = getOrderedArgument(args);
	if (commandsFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# commandsFileName = "+commandsFileName);

	commands = new CommandList(commandsFileName);
	
	modelFileName = getOrderedArgument(args);
	if (modelFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# initial model file name is "+modelFileName);

	String seedString = getOrderedArgument(args);
	if (seedString== null) throw new RuntimeException(errorMessage);
	int seed = (new Integer(seedString)).intValue();
	System.out.println("# seed is "+seed);
	initRandom(seed);
    }
}
