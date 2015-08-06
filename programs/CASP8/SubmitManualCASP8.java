package programs.CASP8;
import meshi.molecularElements.Protein;
import meshi.util.MeshiProgram;


public class SubmitManualCASP8 extends MeshiProgram { 

    private static String proteinFileName = null; 
    private static String targetName = null;
    private static String parentID = null;  
    private static String outputFile = null;  
    private static int model = -1;  
 
    public static void main(String[] args) throws Exception {
		init(args);
		Protein prot = new Protein(proteinFileName);
		prot.printCaspFormat(outputFile , targetName , "5472-6109-2149" , model , parentID , "Energy-based ranking of server models. Specialized loop-building on homology modeling cases. Refinement with minimization under the KB01 potential (Summa & Levitt, 2007)");
	}



    /** ================================= init =========================================
     *
     *A static function for parsing of the command line arguments and assigning the 
     *variables commandsFileName, modelFileName and randomNumberSeed with the right inputs. Note that this
     *static method is using parsing functions such as getOrderedArguments that are defined in MeshiProgram
     *that MinimizeProtein inherits.
     **/
     
    protected static void init(String[] args) {

	String errorMessage = ("\n                  ******************\n"+
			       "Usage java SubmitManualCASP8 <protein file name> <target#> <model#> <parent pdb ID> <output file name> \n"+
			       "                    ******************\n");
			      
	if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));
	proteinFileName = getOrderedArgument(args);
	if (proteinFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# proteinFileName = "+proteinFileName);

	targetName = getOrderedArgument(args);
	if (targetName == null) throw new RuntimeException(errorMessage);
	System.out.println("# targetName = "+targetName);

	String modelString = getOrderedArgument(args);
	if (modelString == null) throw new RuntimeException(errorMessage);
	model = (new Integer(modelString)).intValue();
	System.out.println("# model = "+model);

	parentID = getOrderedArgument(args);
	if (parentID == null) throw new RuntimeException(errorMessage);
	System.out.println("# Parent ID = "+parentID);

	outputFile = getOrderedArgument(args);
	if (outputFile == null) throw new RuntimeException(errorMessage);
	System.out.println("# OutputFileName = "+outputFile);
	
	initRandom(0);
    }

}
