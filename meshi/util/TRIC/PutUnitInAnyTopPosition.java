package meshi.util.TRIC;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Vector;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.overlap.Overlap;

public class PutUnitInAnyTopPosition extends MeshiProgram implements Residues { 
   
    private static String unitLetter = null;
 
    private static int positionNumber = -1;  

    private static int firstResOfCut = -1;  

    private static int lastResOfCut = -1;  

    private static String outputFile = null;  
 
    public static void main(String[] args) {
	init(args); 
	
	try{
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
		Protein stat = new ExtendedAtomsProtein("pos_"+positionNumber+".pdb",DO_NOT_ADD_ATOMS);
		Protein move = new ExtendedAtomsProtein("Atemp_HM_"+unitLetter.charAt(0)+".pdb",DO_NOT_ADD_ATOMS);
		alignByEquatorialDomains(stat.atoms(), 'I', move.atoms(), unitLetter.charAt(0));
		for (int atomC=0 ; atomC<move.atoms().size() ; atomC++) {
			move.atoms().atomAt(atomC).setChain(unitLetter);
		}
		for (int atomC=0 ; atomC<move.atoms().size() ; atomC++) {
			if ((move.atoms().atomAt(atomC).residueNumber()<firstResOfCut) ||
					(move.atoms().atomAt(atomC).residueNumber()>lastResOfCut)) {
				bw.write(move.atoms().atomAt(atomC) + "\n");
			}
		}
		bw.write("TER\n");
		bw.write("END\n");
		bw.close();
	}
	catch(Exception e) {
	    throw new RuntimeException(e.getMessage());
	}    	
}


    /**
     * This procedure will align in space the 'move' protein relative to the equatorial domain
     * of the 'stat' protein. You need to give the Chain IDs so that the residues to be aligned
     * are figured out. Themosome A unit is I, and Thermosome B unit is J, and the 1Q3R is K. 
     * The units of the of Tric are just A,B,G,D,E,H,Q,Z. 
     */
    public static void alignByEquatorialDomains(AtomList statProt, char statProtChainID,
    		AtomList moveProt, char moveProtChainID) {
    	// The ranges:
    	int[][] Q3R = {{25,35},
    			{98,114},
    			{124,140}}; 
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
    	int[][] TCPE = {{34,44},
    			{107,123},
    			{133,149}};
    	
    	int[][] statRange;
    	int[][] moveRange;

    	switch (moveProtChainID) {
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
    	case 'K': 
    		moveRange = Q3R;
    		break;
    	default:
    		throw new RuntimeException("Invalid unit letter {A,B,G,D,H,Q,Z,E,I,J}");
    	}

    	switch (statProtChainID) {
    	case 'A': 
    		statRange = TCPA;
    		break;
    	case 'B': 
    		statRange = TCPB;
    		break;
    	case 'G': 
    		statRange = TCPG;
    		break;
    	case 'D': 
    		statRange = TCPD;
    		break;
    	case 'H': 
    		statRange = TCPH;
    		break;
    	case 'Q': 
    		statRange = TCPQ;
    		break;
    	case 'Z': 
    		statRange = TCPZ;
    		break;
    	case 'E': 
    		statRange = TCPE;
    		break;
    	case 'I': 
    		statRange = ThermoA;
    		break;
    	case 'J': 
    		statRange = ThermoB;
    		break;
    	case 'K': 
    		statRange = Q3R;
    		break;
    	default:
    		throw new RuntimeException("Invalid unit letter {A,B,G,D,H,Q,Z,E,I,J,K}");
    	}

    	
    	Vector<double[]> statArray = new Vector<double[]>();
    	Vector<double[]> moveArray = new Vector<double[]>();
    	int countingRefAtoms = 0;
    	for (int rangeC=0 ; rangeC<moveRange.length ; rangeC++) {
    		for (int resC=moveRange[rangeC][0] ; resC<=moveRange[rangeC][1] ; resC++) {
    			double[] tmpXYZ_move = new double[3];
    			Atom atom_move = moveProt.findAtomInList("CA", resC);
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
    			Atom atom_stat = statProt.findAtomInList("CA", resC);
    			tmpXYZ_stat[0] =atom_stat.x();
    			tmpXYZ_stat[1] =atom_stat.y();
    			tmpXYZ_stat[2] =atom_stat.z();
    			statArray.add(tmpXYZ_stat);
    		}
    	}

    	for (int atomC=0 ; atomC<moveProt.size() ; atomC++) {
    		double[] tmpXYZ_move = new double[3];
    		Atom atom_move = moveProt.atomAt(atomC);
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
    		throw new RuntimeException("The rms value is too large (should be less than 0.6)!");
    	}

    	
    	for (int atomC=0 ; atomC<moveProt.size() ; atomC++) {
    		moveProt.atomAt(atomC).setXYZ(moveAr[0][countingRefAtoms+atomC], 
    				moveAr[1][countingRefAtoms+atomC], moveAr[2][countingRefAtoms+atomC]);
    	}    	
    }
    
    /**
     * 
     * @param chainID - Must be {A,B,G,D,E,H,Q,Z}
     * @param pos - -8,-7,...,-1,1,2,...8
     * 
     * returns an AtomList of the unit in that position.
     * 
     **/
    public static AtomList putAnyUnitInAnyPosition(String chainID, int pos) {
    	if (pos>0) { // upper ring
			AtomList stat = new AtomList("pos_"+pos+".pdb");
			AtomList move = null;
			move = new AtomList("HM_"+chainID.charAt(0)+".pdb");
			alignByEquatorialDomains(stat, 'K', move, chainID.charAt(0));
			return move;
    	}
    	else { // lower ring
			AtomList stat = new AtomList("pos_M"+Math.abs(pos)+".pdb");
			AtomList move = null;
			move = new AtomList("HM_"+chainID.charAt(0)+".pdb");
			alignByEquatorialDomains(stat, 'K', move, chainID.charAt(0));
			return move;
    	}    	
    }
    

    public static AtomList buildFullComplex(String upperSeq , String lowerSeq) {
		AtomList fullComplex = new AtomList();
		// Doing the top ring
		for (int unit=0 ; unit<8 ; unit++) {
			Atom.resetNumberOfAtoms();
			AtomList move = putAnyUnitInAnyPosition(upperSeq.charAt(unit)+"", (unit+1));
			move.setChain(""+upperSeq.charAt(unit));
			fullComplex.add(move);					
		}
		// Doing the bottom ring
		for (int unit=0 ; unit<8 ; unit++) {
			Atom.resetNumberOfAtoms();
			AtomList move = putAnyUnitInAnyPosition(lowerSeq.charAt(unit)+"", -(unit+1));
			move.setChain(""+matchChainLetterOnBottomRing(lowerSeq.charAt(unit)));
			fullComplex.add(move);					
		}
		return fullComplex;    		
	}
    
    
    public static char matchChainLetterOnBottomRing(char upperChain) {
    	char outputChain;
    	switch (upperChain) {
    	case 'A': 
    		outputChain = 'I';
    		break;
    	case 'B': 
    		outputChain = 'J';
    		break;
    	case 'G': 
    		outputChain = 'K';
    		break;
    	case 'D': 
    		outputChain = 'L';
    		break;
    	case 'E': 
    		outputChain = 'M';
    		break;
    	case 'H': 
    		outputChain = 'N';
    		break;
    	case 'Q': 
    		outputChain = 'O';
    		break;
    	case 'Z': 
    		outputChain = 'P';
    		break;
    	default:
    		throw new RuntimeException("Invalid unit letter {A,B,G,D,H,Q,Z,E}");
    	}
    	return outputChain;
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
			       "Usage java  PutUnitInAnyTopPosition <Unit letter: AZBGQDEH> <position [1-8]> <cutting from> <cutting to> <output file for the complex> \n"+
			       "                    ******************\n");
			      
	unitLetter = getOrderedArgument(args);
	if (unitLetter == null) throw new RuntimeException(errorMessage);
	unitLetter = unitLetter.trim();
	System.out.println("# Moving unit: " + unitLetter);

	String positionString = getOrderedArgument(args);
	if (positionString == null)  throw new RuntimeException(errorMessage);
	positionNumber = Integer.parseInt(positionString);
	System.out.println("# To position: " + positionNumber);

	String tmpStr = getOrderedArgument(args);
	if (tmpStr == null) throw new RuntimeException(errorMessage);
	firstResOfCut = Integer.valueOf(tmpStr).intValue();
	System.out.println("# Cutting from: " + firstResOfCut);
	
	tmpStr = getOrderedArgument(args);
	if (tmpStr == null) throw new RuntimeException(errorMessage);
	lastResOfCut = Integer.valueOf(tmpStr).intValue();
	System.out.println("# Cutting to: " + lastResOfCut);

	outputFile = getOrderedArgument(args);
	if (outputFile == null) throw new RuntimeException(errorMessage);
	System.out.println("# Output file name is "+outputFile);
	
	initRandom(999);
    }
}
