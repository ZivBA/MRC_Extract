package programs;

import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.molecularElements.Protein;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;


public class Tric_put_2_proteins_in_a_file extends MeshiProgram implements Residues { 
   
    private static String leftProtFile = null;
 
    private static String rightProtFile = null;  

    private static String outputFile = null;  
 
    public static void main(String[] args) {
	init(args); 


	Protein left = new ExtendedAtomsProtein(leftProtFile,DO_NOT_ADD_ATOMS);
	Protein right = new ExtendedAtomsProtein(rightProtFile,DO_NOT_ADD_ATOMS);
	
	for (int cc=0 ; cc<left.atoms().size() ; cc++)
		left.atoms().atomAt(cc).setChain("A");
	for (int cc=0 ; cc<right.atoms().size() ; cc++)
		right.atoms().atomAt(cc).setChain("B");

	// Outputing
	try{
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
		for (int atomC=0 ; atomC<left.atoms().size() ; atomC++) {
			bw.write(left.atoms().atomAt(atomC) + "\n");
		}
		bw.write("TER\n");
		for (int atomC=0 ; atomC<right.atoms().size() ; atomC++) {
			bw.write(right.atoms().atomAt(atomC) + "\n");
		}
		bw.write("TER\nEND\n");
	    bw.close();
	}
	catch(Exception e) {
	    throw new RuntimeException(e.getMessage());
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

	String errorMessage = ("\n                  ******************\n"+
			       "Usage java  Tric_put_2_proteins_in_a_file <lesft prot> <right prot> <output file name> \n"+
			       "                    ******************\n");
			      
	leftProtFile = getOrderedArgument(args);
	if (leftProtFile == null) throw new RuntimeException(errorMessage);
	System.out.println("# Left protein is: "+leftProtFile);

	rightProtFile = getOrderedArgument(args);
	if (rightProtFile == null) throw new RuntimeException(errorMessage);
	System.out.println("# Right protein is: "+rightProtFile);

	outputFile = getOrderedArgument(args);
	if (outputFile == null) throw new RuntimeException(errorMessage);
	System.out.println("# Output file name is "+outputFile);
	
	initRandom(999);
    }
}
