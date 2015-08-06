package programs;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Vector;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.molecularElements.Atom;
import meshi.molecularElements.Protein;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.overlap.Overlap;


public class Tric_Align_Equator_Domains extends MeshiProgram implements Residues { 
   
    private static String protMoveFile = null;
 
    private static String protStatFile = null;  

    private static char unitLetter = 'X';  

    private static String outputFile = null;  
 
    public static void main(String[] args) {
	init(args); 


	Protein stat = new ExtendedAtomsProtein(protStatFile,DO_NOT_ADD_ATOMS);
	Protein move = new ExtendedAtomsProtein(protMoveFile,DO_NOT_ADD_ATOMS);
	
	// The ranges:
	int[][] ThermoA = {{24,34},
			{97,113},
			{123,139}}; 
	int[][] ThermoB = {{23,33},
			{96,112},
			{122,138}}; 
	int[][] TCPA = {{18,28},
			{91,107},
			{117,133}};
	int[][] TCPB = {{25,35},
			{100,116},
			{126,142}};
	int[][] TCPG = {{23,33},
			{96,112},
			{122,138}};
	int[][] TCPD = {{37,47},
			{110,126},
			{136,152}};
	int[][] TCPH = {{22,32},
			{95,111},
			{121,137}};
	int[][] TCPQ = {{29,39},
			{102,118},
			{128,144}};
	int[][] TCPZ = {{20,30},
			{93,109},
			{119,135}};
	int[][] TCPE = {{20,30},
			{93,109},
			{119,135}};
	
	int[][] statRange = ThermoA;
	int[][] moveRange;
	switch (unitLetter) {
	case 'A': 
		moveRange = TCPA;
		break;
	case 'B': 
		moveRange = TCPB;
		break;
	case 'G': 
		moveRange = TCPG;
		break;
	case 'D': 
		moveRange = TCPD;
		break;
	case 'H': 
		moveRange = TCPH;
		break;
	case 'Q': 
		moveRange = TCPQ;
		break;
	case 'Z': 
		moveRange = TCPZ;
		break;
	case 'E': 
		moveRange = TCPE;
		break;
	case 'I': 
		moveRange = ThermoA;
		break;
	case 'J': 
		moveRange = ThermoB;
		break;
	default:
		throw new RuntimeException("Invalid unit letter {A,B,G,D,H,Q,Z,E,I,J}");
	}


	Vector<double[]> statArray = new Vector<double[]>();
	Vector<double[]> moveArray = new Vector<double[]>();
	int countingRefAtoms = 0;
	for (int rangeC=0 ; rangeC<moveRange.length ; rangeC++) {
		for (int resC=moveRange[rangeC][0] ; resC<=moveRange[rangeC][1] ; resC++) {
			double[] tmpXYZ_move = new double[3];
			Atom atom_move = move.residue(resC).ca();
			tmpXYZ_move[0] =atom_move.x();
			tmpXYZ_move[1] =atom_move.y();
			tmpXYZ_move[2] =atom_move.z();
			moveArray.add(tmpXYZ_move);
			countingRefAtoms++;
		}
	}
	for (int rangeC=0 ; rangeC<statRange.length ; rangeC++) {
		for (int resC=statRange[rangeC][0] ; resC<=statRange[rangeC][1] ; resC++) {
			double[] tmpXYZ_stat = new double[3];
			Atom atom_stat = stat.residue(resC).ca();
			tmpXYZ_stat[0] =atom_stat.x();
			tmpXYZ_stat[1] =atom_stat.y();
			tmpXYZ_stat[2] =atom_stat.z();
			statArray.add(tmpXYZ_stat);
		}
	}

	for (int atomC=0 ; atomC<move.atoms().size() ; atomC++) {
		double[] tmpXYZ_move = new double[3];
		Atom atom_move = move.atoms().atomAt(atomC);
		tmpXYZ_move[0] =atom_move.x();
		tmpXYZ_move[1] =atom_move.y();
		tmpXYZ_move[2] =atom_move.z();
		moveArray.add(tmpXYZ_move);
	}
	
	double[][] statAr = new double[3][statArray.size()];
	double[][] moveAr = new double[3][moveArray.size()];
	int[] common = new int[countingRefAtoms];
	for (int c=0 ; c<statAr[0].length ; c++) {
		statAr[0][c] = statArray.get(c)[0];
		statAr[1][c] = statArray.get(c)[1];
		statAr[2][c] = statArray.get(c)[2];
	}
	for (int c=0 ; c<moveAr[0].length ; c++) {
		moveAr[0][c] = moveArray.get(c)[0];
		moveAr[1][c] = moveArray.get(c)[1];
		moveAr[2][c] = moveArray.get(c)[2];
	}
	for (int c=0 ; c<common.length ; c++) {
		common[c] = c;
	}
	
	
	// Calculating the overlap
	double rmsVal = Overlap.rmsPartialAltRMS(statAr, moveAr, common);
	if (Math.abs(rmsVal)>0.6) {
		throw new RuntimeException("The rms value is not 0.0!");
	}
	else {
		System.out.println(rmsVal);
	}

	
	for (int atomC=0 ; atomC<move.atoms().size() ; atomC++) {
		move.atoms().atomAt(atomC).setXYZ(moveAr[0][countingRefAtoms+atomC], 
				moveAr[1][countingRefAtoms+atomC], moveAr[2][countingRefAtoms+atomC]);
	}
	

	// Outputing
	try{
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
		for (int atomC=0 ; atomC<move.atoms().size() ; atomC++) {
			bw.write(move.atoms().atomAt(atomC) + "\n");
		}
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
			       "Usage java  Tric_Align_Equator_Domains <stat prot> <moving prot> <letter> <output file name> \n"+
			       "                    ******************\n");
			      
	protStatFile = getOrderedArgument(args);
	if (protStatFile == null) throw new RuntimeException(errorMessage);
	System.out.println("# The stationary protein is: "+protStatFile);

	protMoveFile = getOrderedArgument(args);
	if (protMoveFile == null) throw new RuntimeException(errorMessage);
	System.out.println("# The moving protein is: "+protMoveFile);

	String tmp = getOrderedArgument(args);
	if (tmp == null) throw new RuntimeException(errorMessage);
	unitLetter = tmp.charAt(0);
	System.out.println("# Unit letter is: "+unitLetter);

	outputFile = getOrderedArgument(args);
	if (outputFile == null) throw new RuntimeException(errorMessage);
	System.out.println("# Output file name is "+outputFile);

	
	initRandom(999);
    }
}
