package programs;

import meshi.applications.prediction.GDTcalculator;
import meshi.molecularElements.AtomList;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.KeyWords;
import meshi.util.MeshiProgram;


public class CalcRMS extends MeshiProgram implements Residues, AtomTypes , KeyWords { 

	private static String modelFileName = null;
	private static String referenceFileName = null;  

	public static void main(String[] args) {
		init(args); 
		AtomList fullReference = (new AtomList(referenceFileName)).CAFilter();
		AtomList fullProtein = (new AtomList(modelFileName)).CAFilter();
		System.out.println(fullReference.rms(fullProtein));
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
				"java CalcRMS <reference prot> <model prot> \n\n"+
		"                    ******************\n");

		initRandom(0);

		referenceFileName = getOrderedArgument(args);
		if (referenceFileName == null) throw new RuntimeException(errorMessage);

		modelFileName = getOrderedArgument(args);
		if (modelFileName == null) throw new RuntimeException(errorMessage);
	} // of init
} // Of CalcGDT
