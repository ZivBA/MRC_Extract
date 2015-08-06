package programs;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.StringTokenizer;

import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;


public class Tric_Extract_Pairwise_Alignment extends MeshiProgram { 
   
    private static String allAlignments = null;
 
    private static String prefixTop = null;  

    private static String prefixBottom = null;  

    private static String startTop = null;  

    private static String startBottom = null;  

    private static String output = null;  
 
    public static void main(String[] args) {
	init(args); 

    String topString = "";
    String bottomString = "";
    
    String[] clustal = File2StringArray.f2a(allAlignments);
    for (int c=0 ; c<clustal.length ; c++) {
    	StringTokenizer st = new StringTokenizer(clustal[c]);
    	if (st.countTokens()>=2) {
    		String prefix = st.nextToken();
    		String seq = st.nextToken();
    		if (prefix.equals(prefixTop)) { 
    			topString += seq;
    		}
    		if (prefix.equals(prefixBottom)) { 
    			bottomString += seq;
    		}
    	}
    }
    
    // Outputing
	try{
		BufferedWriter bw = new BufferedWriter(new FileWriter(output));
		bw.write("template.pdb\n");
		bw.write(startTop +"\n" + startBottom +"\n");
		bw.write(topString + "\n  \n" + bottomString + "\n");
		bw.write("homology.pdb\n\n");
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
			       "Usage java  Tric_Extract_Pairwise_Alignment <Filename with all alignments> <prefix top> <prefix bottom> <top residue start> <bottom residue start> <output file name> \n"+
			       "                    ******************\n");
			      
	allAlignments = getOrderedArgument(args);
	if (allAlignments == null) throw new RuntimeException(errorMessage);
	System.out.println("# All-Alignments file name is "+allAlignments);

	prefixTop = getOrderedArgument(args);
	if (prefixTop == null) throw new RuntimeException(errorMessage);
	System.out.println("# prefix on top is "+prefixTop);

	prefixBottom = getOrderedArgument(args);
	if (prefixBottom == null) throw new RuntimeException(errorMessage);
	System.out.println("# prefix at bottom is "+prefixBottom);

	startTop = getOrderedArgument(args);
	if (startTop == null) throw new RuntimeException(errorMessage);
	System.out.println("# start on top is "+startTop);

	startBottom = getOrderedArgument(args);
	if (startBottom == null) throw new RuntimeException(errorMessage);
	System.out.println("# start at bottom is "+startBottom);

	output = getOrderedArgument(args);
	if (output == null) throw new RuntimeException(errorMessage);
	System.out.println("# Output file name is "+output);

	
	initRandom(999);
    }
}
