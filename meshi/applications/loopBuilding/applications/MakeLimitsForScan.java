package meshi.applications.loopBuilding.applications;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Vector;

import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.MinimizerException;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;

public class MakeLimitsForScan extends MeshiProgram implements Residues, AtomTypes {

	private static String templateFileName = null;  
	
	private static String outFileName = null;

	private static int loopLength = -1;

	private static int shift = -1;

	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		
		Protein template = new Protein(templateFileName, new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		
		int[][] gaps = getScanLimits(template, loopLength, shift);

		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFileName));
			for (int c=0 ; c<gaps.length ; c++) {
				bw.write(gaps[c][0] + " " + gaps[c][1] + " 999 \n");
				System.out.println("Gap " + c + ": " + gaps[c][0] + " " + gaps[c][1] + " 999 ");							
			}
			bw.close();			
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}				
	}
	
	
	protected static int[][] getScanLimits(Protein template, int loopLength_local, int shift_local) {
		Vector<int[]> resultVec = new Vector<int[]>();
		int startLoop = template.firstResidue();
		int endLoop = startLoop + loopLength_local - 1;
		int lastRes = template.lastResidue();
		int gapCounter = 0;
		while (startLoop<(lastRes-4)) {
			int outStart = startLoop;
			int outEnd = Math.min(endLoop, lastRes);
			if (outEnd!=lastRes-1) { // This will ensure the loop builder wont crash because there is no 2 residues for the terminal overlap.
				int[] tmp = {outStart, outEnd};
				resultVec.add(tmp);
				gapCounter++;
			}
			startLoop += shift_local;
			endLoop += shift_local;
		}
		int[][] outArray = new int[gapCounter][];
		for (int c=0 ; c<gapCounter ; c++) {
			outArray[c] = resultVec.get(c);
		}
		return outArray;		
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
				"Usage java -Xmx600m MakeLimitsForScan <model filename> <output filename> <Loop Length> <Shift>\n"+
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));

		templateFileName = getOrderedArgument(args);
		if (templateFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Template file name is: "+templateFileName);

		outFileName = getOrderedArgument(args);
		if (outFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Output to: "+outFileName);

		String tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		loopLength = (new Integer(tmp.trim())).intValue();
		System.out.println("# Will scan with a loop of this length: " + loopLength);
		
		tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		shift = (new Integer(tmp.trim())).intValue();
		System.out.println("# The shift will be of this length: " + shift);

		initRandom(999);
	}	

} // Of AddMissingResidues
