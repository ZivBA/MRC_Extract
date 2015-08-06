package meshi.util.TRIC;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Vector;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.molecularElements.Atom;
import meshi.molecularElements.Protein;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.overlap.Overlap;

public class HMofFullComplex extends MeshiProgram implements Residues { 
   
    private static String arrangementTopRing = null;
 
    private static String arrangementBottomRing = null;  

    private static int firstResOfCut = -1;  

    private static int lastResOfCut = -1;  

    private static String outputFile = null;  
 
    public static void main(String[] args) {
	init(args); 
	
	System.out.println("\n\nThe inter-ring arrangement:\n-------------------------");
	System.out.println(arrangementTopRing);
	System.out.print(arrangementBottomRing.charAt(0));
	for (int unit=7 ; unit>0 ; unit--) { System.out.print(arrangementBottomRing.charAt(unit)); }
	System.out.println("\n\n");
	
	try{
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));

		// Doing the top ring
		for (int unit=0 ; unit<8 ; unit++) {
			Protein stat = new ExtendedAtomsProtein("pos_"+(unit+1)+".pdb",DO_NOT_ADD_ATOMS);
			Protein move = new ExtendedAtomsProtein("Atemp_HM_"+arrangementTopRing.charAt(unit)+".pdb",DO_NOT_ADD_ATOMS);
			alignByEquatorialDomains(stat, 'I', move, arrangementTopRing.charAt(unit));
			for (int atomC=0 ; atomC<move.atoms().size() ; atomC++) {
				move.atoms().atomAt(atomC).setChain(""+arrangementTopRing.charAt(unit));
			}
			for (int atomC=0 ; atomC<move.atoms().size() ; atomC++) {
				if ((move.atoms().atomAt(atomC).residueNumber()<firstResOfCut) ||
						(move.atoms().atomAt(atomC).residueNumber()>lastResOfCut)) {
					bw.write(move.atoms().atomAt(atomC) + "\n");
				}
			}
			bw.write("TER\n");
			Atom.resetNumberOfAtoms();
		}

		// Doing the bottom ring
		for (int unit=0 ; unit<8 ; unit++) {
			Protein stat = new ExtendedAtomsProtein("pos_M"+(unit+1)+".pdb",DO_NOT_ADD_ATOMS);
			Protein move = new ExtendedAtomsProtein("Atemp_HM_"+arrangementBottomRing.charAt(unit)+".pdb",DO_NOT_ADD_ATOMS);
			alignByEquatorialDomains(stat, 'I', move, arrangementBottomRing.charAt(unit));
			// Finding the output chain letter for the lower ring
	    	char outputChain = PutUnitInAnyTopPosition.matchChainLetterOnBottomRing(arrangementBottomRing.charAt(unit));
			for (int atomC=0 ; atomC<move.atoms().size() ; atomC++) {
				move.atoms().atomAt(atomC).setChain(""+outputChain);
			}
			for (int atomC=0 ; atomC<move.atoms().size() ; atomC++) {
				if ((move.atoms().atomAt(atomC).residueNumber()<firstResOfCut) ||
						(move.atoms().atomAt(atomC).residueNumber()>lastResOfCut)) {
					bw.write(move.atoms().atomAt(atomC) + "\n");
				}
			}
			bw.write("TER\n");
			Atom.resetNumberOfAtoms();
		}
		
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
     * are figured out. Themosome A unit is I, and Thermosome B unit is J. The units of the of Tric
     * are just A,B,G,D,E,H,Q,Z. 
     */
    protected static void alignByEquatorialDomains(Protein statProt, char statProtChainID,
    		 Protein moveProt, char moveProtChainID) {
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
    	default:
    		throw new RuntimeException("Invalid unit letter {A,B,G,D,H,Q,Z,E,I,J}");
    	}

    	
    	Vector<double[]> statArray = new Vector<double[]>();
    	Vector<double[]> moveArray = new Vector<double[]>();
    	int countingRefAtoms = 0;
    	for (int rangeC=0 ; rangeC<moveRange.length ; rangeC++) {
    		for (int resC=moveRange[rangeC][0] ; resC<=moveRange[rangeC][1] ; resC++) {
    			double[] tmpXYZ_move = new double[3];
    			Atom atom_move = moveProt.residue(resC).ca();
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
    			Atom atom_stat = statProt.residue(resC).ca();
    			tmpXYZ_stat[0] =atom_stat.x();
    			tmpXYZ_stat[1] =atom_stat.y();
    			tmpXYZ_stat[2] =atom_stat.z();
    			statArray.add(tmpXYZ_stat);
    		}
    	}

    	for (int atomC=0 ; atomC<moveProt.atoms().size() ; atomC++) {
    		double[] tmpXYZ_move = new double[3];
    		Atom atom_move = moveProt.atoms().atomAt(atomC);
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

    	
    	for (int atomC=0 ; atomC<moveProt.atoms().size() ; atomC++) {
    		moveProt.atoms().atomAt(atomC).setXYZ(moveAr[0][countingRefAtoms+atomC], 
    				moveAr[1][countingRefAtoms+atomC], moveAr[2][countingRefAtoms+atomC]);
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
			       "Usage java  HMofFullComplex <Top ring arrangemen like: AZBGQDEH> <Bottom ring arrangemen like: AZBGQDEH> <first Res of cut> <last Res of cut> <output file for the complex> \n"+
			       "                    ******************\n");
			      
	arrangementTopRing = getOrderedArgument(args);
	if ((arrangementTopRing == null) || (arrangementTopRing.length()!=8)) throw new RuntimeException(errorMessage);
	System.out.println("# The top ring arrangement is: " + arrangementTopRing);

	arrangementBottomRing = getOrderedArgument(args);
	if ((arrangementBottomRing == null) || (arrangementBottomRing.length()!=8)) throw new RuntimeException(errorMessage);
	System.out.println("# The bottom ring arrangement is: " + arrangementBottomRing);

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
