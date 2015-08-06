package meshi.util.TRIC;

import java.util.Vector;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.overlap.Overlap;

public class PutUnitInAnyTopPositionYeast extends MeshiProgram implements Residues { 
   

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
    	int[][] ThermoA = null; 
    	int[][] ThermoB = null; 
    	int[][] TCPA = {{26,36}};
    	int[][] TCPB = {{19,29}};
    	int[][] TCPG = {{21,31}};
    	int[][] TCPD = {{21,31}};
    	int[][] TCPE = {{47,57}};
    	int[][] TCPH = {{26,36}};
    	int[][] TCPQ = {{28,38}};
    	int[][] TCPZ = {{19,29}};
    	
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
			move = new AtomList("HM_yeast_"+chainID.charAt(0)+".pdb");
			alignByEquatorialDomains(stat, 'K', move, chainID.charAt(0));
			return move;
    	}
    	else { // lower ring
			AtomList stat = new AtomList("pos_M"+Math.abs(pos)+".pdb");
			AtomList move = null;
			move = new AtomList("HM_yeast_"+chainID.charAt(0)+".pdb");
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


}
