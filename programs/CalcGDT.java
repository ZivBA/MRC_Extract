package programs;

import java.util.Vector;

import meshi.applications.prediction.GDTcalculator;
import meshi.molecularElements.AtomList;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.KeyWords;
import meshi.util.MeshiProgram;


public class CalcGDT extends MeshiProgram implements Residues, AtomTypes , KeyWords { 

    private static String modelFileName = null;
    private static String referenceFileName = null;  
    private static double cutoff1 = 0.5;
    private static double cutoff2 = 1.0;
    private static double cutoff3 = 2.0;
    private static double cutoff4 = 4.0;
    private static Vector<int[]> ranges;
    
 
    public static void main(String[] args) {
	init(args); 
	AtomList fullReference = (new AtomList(referenceFileName)).CAFilter();
	AtomList fullProtein = (new AtomList(modelFileName)).CAFilter();
	AtomList ref = new AtomList();
	AtomList model = new AtomList();
	if (ranges.isEmpty()) { // all the protein
		for (int c=0 ; c<fullReference.size() ; c++) {
			ref.add(fullReference.atomAt(c));
		}
		for (int c=0 ; c<fullProtein.size() ; c++) {
			model.add(fullProtein.atomAt(c));
		}
	}
	else {  // ranges
		for (int rangeCounter=0 ; rangeCounter<ranges.size() ; rangeCounter++) {
			for (int c=0 ; c<fullReference.size() ; c++) {
				if ((fullReference.atomAt(c).residueNumber()>=ranges.get(rangeCounter)[0]) &&
						(fullReference.atomAt(c).residueNumber()<=ranges.get(rangeCounter)[1])) {
					ref.add(fullReference.atomAt(c));
				}
			}
			for (int c=0 ; c<fullProtein.size() ; c++) {
				if ((fullProtein.atomAt(c).residueNumber()>=ranges.get(rangeCounter)[0]) &&
						(fullProtein.atomAt(c).residueNumber()<=ranges.get(rangeCounter)[1])) {
					model.add(fullProtein.atomAt(c));
				}
			}
		}
	}
	if (ref.size() == model.size()) {
		System.out.println("GDT: " + GDTcalculator.gdt(ref,model,cutoff1, cutoff2, cutoff3, cutoff4) +	"      RMS: " + ref.getRms(model));
	}
	else {
		System.out.println("GDT: " + GDTcalculator.gdt(ref,model,cutoff1, cutoff2, cutoff3, cutoff4) +	"      RMS: Not computable. The number of CAs in both PDB file is not the same.");
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


	String errorMessage = ("\n                  ******************\n\nUSAGE:\n------\n"+
			"For calculating on all the protein resdiues with default (high accuracy) cutoffs of 0.5 , 1.0 , 2.0 , 4.0:\n" +
			"java -Xmx300m CalcGDT <reference prot> <model prot> \n\n"+
			"For calculating on all the protein resdiues with your cutoffs:\n" +
			"java -Xmx300m CalcGDT <reference prot> <model prot> -c <cutoff 1>  <cutoff 2>  <cutoff 3>  <cutoff 4>\n\n"+
			"For calculating on a specific (not necessary continous) set of the protein resdiues with your cutoffs:\n" +
			"java -Xmx300m CalcGDT <reference prot> <model prot> -c <cutoff 1>  <cutoff 2>  <cutoff 3>  <cutoff 4> -r <from res (fragment 1)> <to res (fragment 1)> <from res (fragment 2)> <to res (fragment 2)> ...\n\n"+
			"                    ******************\n");
			      
	if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));

	initRandom(0);
	
	referenceFileName = getOrderedArgument(args);
	if (referenceFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# reference model file name is "+referenceFileName);
	
	modelFileName = getOrderedArgument(args);
	if (modelFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# initial model file name is "+modelFileName);
	
	ranges = new Vector<int[]>();
	
	if (!getFlag("-c",args))
		return;
	else {
		String tmp = "";
		
		tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		cutoff1 = (new Double(tmp)).doubleValue();
		System.out.println("# Cutoff 1 is " + cutoff1);
	
		tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		cutoff2 = (new Double(tmp)).doubleValue();
		System.out.println("# Cutoff 2 is " + cutoff2);
	
		tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		cutoff3 = (new Double(tmp)).doubleValue();
		System.out.println("# Cutoff 3 is " + cutoff3);
	
		tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		cutoff4 = (new Double(tmp)).doubleValue();
		System.out.println("# Cutoff 4 is " + cutoff4);
		
		if (!getFlag("-r",args))
			return;
		else {
			tmp = getOrderedArgument(args);
			
			while (tmp != null) {
				int[] tmpRange = new int[2];
				tmpRange[0] = (new Double(tmp)).intValue();
				tmp = getOrderedArgument(args);
				if (tmp==null)
					throw new RuntimeException("Error reading the command line. The number of residue range limiters is not even");
				tmpRange[1] = (new Integer(tmp)).intValue();
				if (tmpRange[1]<tmpRange[0])  
					throw new RuntimeException("the low limit of a range: " + tmpRange[0] + " is higher than the high limit: " + tmpRange[1]);
				ranges.add(tmpRange);
				System.out.println("# Restricting calculations also to this range: ["+ tmpRange[0] +"," + tmpRange[1]+"]");
				tmp = getOrderedArgument(args);
			}
		}
    }
} // of init
} // Of CalcGDT
