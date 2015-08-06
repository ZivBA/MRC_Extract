package meshi.applications.TriC;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.StringTokenizer;

import meshi.util.file.File2StringArray;

public class TricAlignment {
	
	private String Q3R = "";
	private String thermoA = "";
	private String thermoB = "";
	private String A = "";
	private String B = "";
	private String G = "";
	private String D = "";
	private String E = "";
	private String H = "";
	private String Q = "";
	private String Z = "";
	private String W = "";

	public TricAlignment() {
		String[] clustal = File2StringArray.f2a("C:\\Users\\Nir\\TRiC\\Organisms\\TRiC_11_2_2011_alignment_1line");
		for (int c=0 ; c<clustal.length ; c++) {
			StringTokenizer st = new StringTokenizer(clustal[c]);
			if (st.countTokens()>=2) {
				String prefix = st.nextToken();
				String seq = st.nextToken();
				if (prefix.equals("1Q3R=A")) { 
					Q3R += seq;
				}
				if (prefix.equals("1A6D=A")) { 
					thermoA += seq;
				}
				if (prefix.equals("1A6D=B")) { 
					thermoB += seq;
				}
				if (prefix.equals("TRiC=A")) { 
					A += seq;
				}
				if (prefix.equals("TRiC=B")) { 
					B += seq;
				}
				if (prefix.equals("TRiC=G")) { 
					G += seq;
				}
				if (prefix.equals("TRiC=D")) { 
					D += seq;
				}
				if (prefix.equals("TRiC=E")) { 
					E += seq;
				}
				if (prefix.equals("TRiC=H")) { 
					H += seq;
				}
				if (prefix.equals("TRiC=Q")) { 
					Q += seq;
				}
				if (prefix.equals("TRiC=Z")) { 
					Z += seq;
				}
				if (prefix.equals("TRiC=W")) { 
					W += seq;
				}
			}
		}
	}
	
	/*
	 * Position in the alignment (0-based)
	 */
	private int getPositionInSeq(char refUnit, int res) {
    	String refString = null;
		switch (refUnit) {
    	case 'A': 
    		refString = A;
    		break;
    	case 'B': 
    		refString = B;
    		break;
    	case 'G': 
    		refString = G;
    		break;
    	case 'D': 
    		refString = D;
    		break;
    	case 'E': 
    		refString = E;
    		break;
    	case 'H': 
    		refString = H;
    		break;
    	case 'Q': 
    		refString = Q;
    		break;
    	case 'Z': 
    		refString = Z;
    		break;
    	case 'W': 
    		refString = W;
    		break;
    	case 'I': 
    		refString = thermoA;
    		break;
    	case 'J': 
    		refString = thermoB;
    		break;
    	case 'K': 
    		refString = Q3R;
    		break;
    	default:
    		throw new RuntimeException("Invalid unit letter {A,B,G,D,H,Q,Z,W,E,I,J,K}");
    	}
		int posInSeq = 0;
		int resCounter = 0;
		for ( ; resCounter<res ; posInSeq++) {
			if (refString.charAt(posInSeq)!='-')
				resCounter++;
		}
		
		return posInSeq-1;		
	}
	
	/*
	 * translate residues between units.
	 */
	public int getNewResNum(char refUnit, int res, char queryUnit) {
		int posInSeq = getPositionInSeq(refUnit, res); 
  
		String qString = null;
		switch (queryUnit) {
    	case 'A': 
    		qString = A;
    		break;
    	case 'B': 
    		qString = B;
    		break;
    	case 'G': 
    		qString = G;
    		break;
    	case 'D': 
    		qString = D;
    		break;
    	case 'E': 
    		qString = E;
    		break;
    	case 'H': 
    		qString = H;
    		break;
    	case 'Q': 
    		qString = Q;
    		break;
    	case 'Z': 
    		qString = Z;
    		break;
    	case 'W': 
    		qString = W;
    		break;
    	case 'I': 
    		qString = thermoA;
    		break;
    	case 'J': 
    		qString = thermoB;
    		break;
    	case 'K': 
    		qString = Q3R;
    		break;
    	default:
    		throw new RuntimeException("Invalid unit letter {A,B,G,D,H,Q,Z,E,W,I,J,K}");
    	}
		
		int resCounter = 0;
		for (int sCounter=0 ; sCounter<=posInSeq ; sCounter++) {
			if (qString.charAt(sCounter)!='-')
				resCounter++;
		}
		
		return resCounter;		
	}
	
	
	/*
	 * AA letter according to residue
	 */
	public char getAAinSeq(char unit, int res) {
		int posInSeq = getPositionInSeq(unit, res); 
		
		switch (unit) {
    	case 'A': 
    		return A.charAt(posInSeq);
    	case 'B': 
    		return B.charAt(posInSeq);
    	case 'G': 
    		return G.charAt(posInSeq);
    	case 'D': 
    		return D.charAt(posInSeq);
    	case 'E': 
    		return E.charAt(posInSeq);
    	case 'H': 
    		return H.charAt(posInSeq);
    	case 'Q': 
    		return Q.charAt(posInSeq);
    	case 'Z': 
    		return Z.charAt(posInSeq);
    	case 'W': 
    		return W.charAt(posInSeq);
    	case 'I': 
    		return thermoA.charAt(posInSeq);
    	case 'J': 
    		return thermoB.charAt(posInSeq);
    	case 'K': 
    		return Q3R.charAt(posInSeq);
    	default:
    		throw new RuntimeException("Invalid unit letter {A,B,G,D,H,Q,Z,W,E,I,J,K}");
    	}		
	}
	
	/*
	 * translate between units.
	 */
	public char getAAinSeq(char refUnit, int res, char queryUnit) {
		int newResNum = getNewResNum(refUnit, res, queryUnit);
		return getAAinSeq(queryUnit, newResNum);
	}


	public static void main(String[] args) {
		TricAlignment alignments = new TricAlignment();
//		System.out.println(alignments.W);
//		System.out.println();
//		System.out.println(alignments.Q3R);

		System.out.println("1Q3R=A  " + alignments.Q3R);
		System.out.println("1A6D=A  " + alignments.thermoA);
		System.out.println("1A6D=B  " + alignments.thermoB);
		System.out.println("TRiC=A  " + alignments.A);
		System.out.println("TRiC=B  " + alignments.B);
		System.out.println("TRiC=G  " + alignments.G);
		System.out.println("TRiC=D  " + alignments.D);
		System.out.println("TRiC=E  " + alignments.E);
		System.out.println("TRiC=H  " + alignments.H);
		System.out.println("TRiC=Q  " + alignments.Q);
		System.out.println("TRiC=Z  " + alignments.Z);
		System.out.println("TRiC=W  " + alignments.W);
		
		
		
//		try{
//			BufferedWriter bw = new BufferedWriter(new FileWriter(args[0]));
//			bw.write("C:\\Users\\Nir\\TriC\\HM_models_of_units\\1A6D\\Btemplate\\template.pdb\n");
//			bw.write("1\n1\n");
//			String tricString = null;
//			switch (args[1].charAt(0)) {
//	    	case 'A': 
//	    		tricString = alignments.A;
//	    		break;
//	    	case 'B': 
//	    		tricString = alignments.B;
//	    		break;
//	    	case 'G': 
//	    		tricString = alignments.G;
//	    		break;
//	    	case 'D': 
//	    		tricString = alignments.D;
//	    		break;
//	    	case 'E': 
//	    		tricString = alignments.E;
//	    		break;
//	    	case 'H': 
//	    		tricString = alignments.H;
//	    		break;
//	    	case 'Q': 
//	    		tricString = alignments.Q;
//	    		break;
//	    	case 'Z': 
//	    		tricString = alignments.Z;
//	    		break;
//	    	default:
//	    		throw new RuntimeException("Invalid unit letter {A,B,G,D,E,H,Q,Z}");
//	    	}
//			bw.write(tricString + "\n  \n" + alignments.thermoB + "\n");
//			bw.write("homology.pdb\n\n");
//		    bw.close();
//		}
//		catch(Exception e) {
//		    throw new RuntimeException(e.getMessage());
//		}    	
	}

	public String getAlignment(String chainID) {
		switch (chainID.charAt(0)) {
    	case 'A': 
    		return A;
    	case 'B': 
    		return B;
    	case 'G': 
    		return G;
    	case 'D': 
    		return D;
    	case 'E': 
    		return E;
    	case 'H': 
    		return H;
    	case 'Q': 
    		return Q;
    	case 'Z': 
    		return Z;
    	case 'W': 
    		return W;
    	case 'I': 
    		return thermoA;
    	case 'J': 
    		return thermoB;
    	case 'K': 
    		return Q3R;
    	default:
    		throw new RuntimeException("Invalid unit letter {A,B,G,D,H,Q,Z,W,E,I,J,K}");
    	}		
	}	
}
