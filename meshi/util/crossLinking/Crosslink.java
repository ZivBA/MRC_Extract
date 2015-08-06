package meshi.util.crossLinking;

import java.util.StringTokenizer;

import javax.swing.text.Position;



public class Crosslink {
	
	private double score = -1; // some sort of confidence scoring. The higher -> higher confidence.	
	private double scoreRapp = -1; 	
	private double scoreKalis = -1; 	
	private String seq1 = null; // The sequence of the first peptide.
	private String seq2 = null; // The sequence of the second peptide.
	private int xl_pos_in_seq1 = -1; // A zero-based residue number in peptide 1 of the crosslinked residue.
	private int xl_pos_in_seq2 = -1; // A zero-based residue number in peptide 2 of the crosslinked residue.
	private String protName1 = null; // The protein name of the first peptide.
	private String protName2 = null; // The protein name of the second peptide.
	private String msFilename = null;
	private int absPos1 = -1; // The (one-based) residue number in the entire sequence of the subunit  	
	private int absPos2 = -1; // The (one-based) residue number in the entire sequence of the subunit

	/**
	 * Building an object describing a single inter-peptide crosslink. The text line should 
	 * be from the abersold lab excel file.
	 */
	public Crosslink(String textLine, int lineFormat) {
		switch (lineFormat) {
		case 0 : readAebersold(textLine); break;
		case 1 : readKalisman(textLine); break;
		case 2 : readRappsilber(textLine); break;
		case 3 : readPhil(textLine); break;		
		case 4 : readSciencePaper(textLine); break;		
		default: throw new RuntimeException("An unknown XL line format.");
		}
	}
	
	/**
	 * "Hardcoding" the cross-link
	 */
	public Crosslink(double score, String seq1, String seq2, int xl_pos_in_seq1, int xl_pos_in_seq2,
			String protName1, String protName2, String msFilename, int absPos1, int absPos2) {
		this.score = score;
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.xl_pos_in_seq1 = xl_pos_in_seq1;
		this.xl_pos_in_seq2 = xl_pos_in_seq2;
		this.protName1 = protName1;
		this.protName2 = protName2;
		this.msFilename = msFilename;
		this.absPos1 = absPos1;
		this.absPos2 = absPos2;
	}
	
	private void readAebersold(String textLine) {
		StringTokenizer st = new StringTokenizer(textLine);
		// Reading column 1.
		String tmpStr = st.nextToken();
		int separator1 = tmpStr.indexOf('-', 0);
		int separator2 = tmpStr.indexOf('-', separator1+1);
		int separator3 = tmpStr.indexOf('-', separator2+1);
		seq1 = tmpStr.substring(0, separator1);
		seq2 = tmpStr.substring(separator1+1, separator2);
		xl_pos_in_seq1 = Integer.parseInt(tmpStr.substring(separator2+2, separator3)) - 1;
		xl_pos_in_seq2 = Integer.parseInt(tmpStr.substring(separator3+2)) - 1;
		// Reading column 2.
		protName1 = st.nextToken();		
		// Reading column 3.
		protName2 = st.nextToken();		
		// Reading column 4.
		absPos1 = Integer.parseInt(st.nextToken());		
//		if (protName1.equals("D")) { // For Bovine
//			absPos1 += 3;
//		}
		// Reading column 5
		absPos2 = Integer.parseInt(st.nextToken());		
//		if (protName2.equals("D")) {  // For Bovine
//			absPos2 += 3;
//		}		
		st.nextToken();
		st.nextToken();
		st.nextToken();
		st.nextToken();
		// Reading column 6.
		scoreRapp = Double.parseDouble(st.nextToken());
		scoreKalis = scoreRapp;
		//		msFilename = st.nextToken();		
	}

	private void readKalisman(String textLine) {
		StringTokenizer st = new StringTokenizer(textLine);
		// Reading column 1.
		protName1 = st.nextToken();		
		// Reading column 2.
		protName2 = st.nextToken();		
		// Reading column 3.
		absPos1 = Integer.parseInt(st.nextToken());		
		if (protName1.equals("Z") & (absPos1>373)) { // The yeast insertion
			absPos1 -= 56;
		}
		// Reading column 4.
		absPos2 = Integer.parseInt(st.nextToken());		
		if (protName2.equals("Z") & (absPos2>373)) { // The insertion
			absPos2 -= 56;
		}
		// Reading column 5.
		scoreRapp = Double.parseDouble(st.nextToken()); // The high/low/bad score
		score = scoreRapp;
		// Reading column 6.
		scoreKalis = Double.parseDouble(st.nextToken()); // The Smsms score
		// Reading column 7.
		st.nextToken(); // The first in 10 score
		// Reading column 8.
		st.nextToken(); // The intensity
		// Reading column 9.
		st.nextToken(); // The charge
		// Reading column 10.
		st.nextToken(); // The Mass
		// Reading column 11.
		String tmpStr = st.nextToken();
		int separator = tmpStr.indexOf('-');
		seq1 = tmpStr.substring(0, separator);
		xl_pos_in_seq1 = Integer.parseInt(tmpStr.substring(separator+1)) - 1;
		// Reading column 12.
		tmpStr = st.nextToken();
		separator = tmpStr.indexOf('-');
		seq2 = tmpStr.substring(0, separator);
		xl_pos_in_seq2 = Integer.parseInt(tmpStr.substring(separator+1)) - 1;
	}

	private void readPhil(String textLine) {
		StringTokenizer st = new StringTokenizer(textLine);
		st.nextToken(); // Reading column 1 - MZ.
		st.nextToken(); // Reading column 2 - Z.
		st.nextToken(); // Reading column 3 - ppm.
		msFilename = st.nextToken(); // Reading column 4 - fraction name.
		st.nextToken(); // Reading column 5 - RT.
		st.nextToken(); // Reading column 6 - S1.
		scoreRapp = Double.parseDouble(st.nextToken()); // Reading column 7 - S2 score_diff.
		st.nextToken(); // Reading column 8 - S3.
		st.nextToken(); // Reading column 9 - S4.
		scoreKalis = Double.parseDouble(st.nextToken()); // Reading column 10 - S5 matched.
		score = -70.0/13.0*scoreRapp + 70.0*23.0/13.0;
		seq1 = st.nextToken(); // Reading column 11 - seq1.
		seq2 = st.nextToken(); // Reading column 12 - seq2.
		st.nextToken(); // Reading column 13 - ??.
		st.nextToken(); // Reading column 14 - ??.
		st.nextToken(); // Reading column 15 - ??.
		st.nextToken(); // Reading column 16 - ??.
		absPos1 = Integer.parseInt(st.nextToken());  // Reading column 17 - residue number in pep1 		
		protName1 = convertPhil2Us(st.nextToken());  // Reading column 18 - prot of pep1 				
		absPos2 = Integer.parseInt(st.nextToken());  // Reading column 17 - residue number in pep2 		
		protName2 = convertPhil2Us(st.nextToken());  // Reading column 18 - prot of pep2 				
		xl_pos_in_seq1 = -99;
		xl_pos_in_seq2 = -99;
	}

	
	private void readRappsilber(String textLine) {
		StringTokenizer st = new StringTokenizer(textLine);
		// Reading column 1.
		protName1 = convertRapp2Us(st.nextToken());		
		// Reading column 2.
		absPos1 = Integer.parseInt(st.nextToken());		
		// Reading column 3.
		protName2 = convertRapp2Us(st.nextToken());		
		// Reading column 4.
		absPos2 = Integer.parseInt(st.nextToken());		
		// Reading column 5.
		String scoreString = st.nextToken(); // The high/low/bad score
		if (scoreString.equals("High")) {
			score = 2;
		}
		else if (scoreString.equals("Low")) {
			score = 1;
		}
		else {
			throw new RuntimeException("This should not happen:  " + scoreString);
		}
		scoreRapp = score;
		seq1 = "";
		xl_pos_in_seq1 = -99;
		seq2 = "";
		xl_pos_in_seq2 = -99;
	}
	

	private void readSciencePaper(String textLine) {
		StringTokenizer st = new StringTokenizer(textLine);
		// Reading column 1.
		protName1 = st.nextToken();		
		// Reading column 2.
		protName2 = st.nextToken();		
		// Reading column 3.
		absPos1 = Integer.parseInt(st.nextToken());		
		// Reading column 4.
		absPos2 = Integer.parseInt(st.nextToken());		
		// Reading column 5.
		scoreRapp = Double.parseDouble(st.nextToken()); // The high/low/bad score
		score = scoreRapp;
		// Reading column 6.
		scoreKalis = Double.parseDouble(st.nextToken()); // The Smsms score
		// Reading column 7.
		st.nextToken(); // Charge
		// Reading column 8.
		st.nextToken(); // PPM mass error
		// Reading column 9.
		st.nextToken(); // Distance in model
		// Reading column 10 - first peptide.
		String tmpStr = st.nextToken();
		int separator = tmpStr.indexOf('-');
		seq1 = tmpStr.substring(0, separator);
		xl_pos_in_seq1 = Integer.parseInt(tmpStr.substring(separator+1)) - 1;
		// Reading column 11 - second peptide.
		tmpStr = st.nextToken();
		separator = tmpStr.indexOf('-');
		seq2 = tmpStr.substring(0, separator);
		xl_pos_in_seq2 = Integer.parseInt(tmpStr.substring(separator+1)) - 1;
		protName1 = convertScience2Us(protName1);   				
		protName2 = convertScience2Us(protName2);   				
		xl_pos_in_seq1 = -99;
		xl_pos_in_seq2 = -99;
	}

	
	private static String convertRapp2Us(String rapp) {
		if (rapp.equals("Rpb1")) {
			return "A";
		}
		if (rapp.equals("Rpb2")) {
			return "B";
		}
		if (rapp.equals("Rpb3")) {
			return "C";
		}
		if (rapp.equals("Rpb4")) {
			return "D";
		}
		if (rapp.equals("Rpb5")) {
			return "E";
		}
		if (rapp.equals("Rpb6")) {
			return "F";
		}
		if (rapp.equals("Rpb7")) {
			return "G";
		}
		if (rapp.equals("Rpb8")) {
			return "H";
		}
		if (rapp.equals("Rpb9")) {
			return "I";
		}
		if (rapp.equals("Rpb10")) {
			return "J";
		}
		if (rapp.equals("Rpb11")) {
			return "K";
		}
		if (rapp.equals("Rpb12")) {
			return "L";
		}
		throw new RuntimeException("Unknown protein name in Rappsilber:  " + rapp);
	}


	private static String convertPhil2Us(String rapp) {
		if (rapp.equals("Rpb1")) {
			return "A";
		}
		if (rapp.equals("Rpb2")) {
			return "B";
		}
		if (rapp.equals("Rpb3")) {
			return "C";
		}
		if (rapp.equals("Rpb4")) {
			return "D";
		}
		if (rapp.equals("Rpb5")) {
			return "E";
		}
		if (rapp.equals("Rpb6")) {
			return "F";
		}
		if (rapp.equals("Rpb7")) {
			return "G";
		}
		if (rapp.equals("Rpb8")) {
			return "H";
		}
		if (rapp.equals("Rpb9")) {
			return "I";
		}
		if (rapp.equals("Rpb10")) {
			return "J";
		}
		if (rapp.equals("Rpb11")) {
			return "K";
		}
		if (rapp.equals("Rpb12")) {
			return "L";
		}
		if (rapp.equals("Med1")) {
			return "M";
		}
		if (rapp.equals("Med2")) {
			return "N";
		}
		if (rapp.equals("Med3")) {
			return "O";
		}
		if (rapp.equals("Med4")) {
			return "P";
		}
		if (rapp.equals("Med5")) {
			return "Q";
		}
		if (rapp.equals("Med6")) {
			return "R";
		}
		if (rapp.equals("Med7")) {
			return "S";
		}
		if (rapp.equals("Med8")) {
			return "T";
		}
		if (rapp.equals("Med9")) {
			return "U";
		}
		if (rapp.equals("Med10")) {
			return "V";
		}
		if (rapp.equals("Med11")) {
			return "W";
		}
		if (rapp.equals("Med14")) {
			return "X";
		}
		if (rapp.equals("Med15")) {
			return "Y";
		}
		if (rapp.equals("Med16")) {
			return "Z";
		}
		if (rapp.equals("Med17")) {
			return "1";
		}
		if (rapp.equals("Med18")) {
			return "2";
		}
		if (rapp.equals("Med19")) {
			return "3";
		}
		if (rapp.equals("Med20")) {
			return "4";
		}
		if (rapp.equals("Med21")) {
			return "5";
		}
		if (rapp.equals("Med22")) {
			return "6";
		}
		throw new RuntimeException("Unknown protein name in Phil\'s data:  " + rapp);
	}

	private static String convertScience2Us(String rapp) {
		if (rapp.equals("Rpb1")) {
			return "A";
		}
		if (rapp.equals("Rpb2")) {
			return "B";
		}
		if (rapp.equals("Rpb3")) {
			return "C";
		}
		if (rapp.equals("Rpb4")) {
			return "D";
		}
		if (rapp.equals("Rpb5")) {
			return "E";
		}
		if (rapp.equals("Rpb6")) {
			return "F";
		}
		if (rapp.equals("Rpb7")) {
			return "G";
		}
		if (rapp.equals("Rpb8")) {
			return "H";
		}
		if (rapp.equals("Rpb9")) {
			return "I";
		}
		if (rapp.equals("Rpb10")) {
			return "J";
		}
		if (rapp.equals("Rpb11")) {
			return "K";
		}
		if (rapp.equals("Rpb12")) {
			return "L";
		}
		if (rapp.equals("Toa1")) {
			return "M";
		}
		if (rapp.equals("Toa2")) {
			return "O";
		}
		if (rapp.equals("TFIIB")) {
			return "P";
		}
		if (rapp.equals("TBP")) {
			return "Q";
		}
		if (rapp.equals("Tfa1")) {
			return "R";
		}
		if (rapp.equals("Tfa2")) {
			return "S";
		}
		if (rapp.equals("Tfg1")) {
			return "U";
		}
		if (rapp.equals("Tfg2")) {
			return "V";
		}
		if (rapp.equals("Tfg3")) {
			return "3";
		}
		if (rapp.equals("Ssl2")) {
			return "1";
		}
		if (rapp.equals("Rad3")) {
			return "Y";
		}
		if (rapp.equals("Tfb1")) {
			return "4";
		}
		if (rapp.equals("Tfb2")) {
			return "W";
		}
		if (rapp.equals("Ssl1")) {
			return "Z";
		}
		if (rapp.equals("Tfb3")) {
			return "5";
		}
		if (rapp.equals("Tfb4")) {
			return "6";
		}
		if (rapp.equals("Ccl1")) {
			return "7";
		}
		if (rapp.equals("Kin28")) {
			return "8";
		}
		if (rapp.equals("Tfb5")) {
			return "X";
		}
		if (rapp.equals("TFIIS")) {
			return "2";
		}
		if (rapp.equals("Sub1")) {
			return "9";
		}
		throw new RuntimeException("Unknown protein name in Science\'s data:  " + rapp);
	}
	
	
	public String toString() {
		return protName1 + " " + absPos1 + " " + seq1 + "  <-XXX->  " +  seq2 + " " + absPos2 + " " + protName2 + "     [" + scoreRapp + "," + scoreKalis + "]\n";
	}
	
	public String toStringKalismanFormat() {
		return protName1() + " " + protName2() + " " + absPos1() + " " + absPos2() + " -1 " + score() + "  -1 -1 -1 -1 " +  seq1() + "-" + (xl_pos_in_seq1()+1) + " " + seq2() + "-" + (xl_pos_in_seq2()+1);
	}

	public String toStringWithFileName() {
		return protName1 + " " + absPos1 + " " + seq1 + "  <-XXX->  " +  seq2 + " " + absPos2 + " " + protName2 + "         " + msFilename + "\n";
	}
	
	public String seq1() {
		return seq1;
	}

	public String seq2() {
		return seq2;
	}

	public int xl_pos_in_seq1() {
		return xl_pos_in_seq1;
	}

	public int xl_pos_in_seq2() {
		return xl_pos_in_seq2;
	}

	public String protName1() {
		return protName1;
	}

	public String protName2() {
		return protName2;
	}

	public void set_MS_file_name(String fileName) {
		msFilename = fileName;
	}

	public String msFilename() {
		return msFilename;
	}

	public int absPos1() {
		return absPos1;
	}

	public int absPos2() {
		return absPos2;
	}

	public double score() {
		return score;
	}

	public double scoreKalis() {
		return scoreKalis;
	}

	public double scoreRapp() {
		return scoreRapp;
	}

	/**
	 * These 2 setters are definitely not good programing. I'm only putting them so that 
	 * I can run the random shuffling experiment, and the 'PuttingXLinPDB'.
	 */
	public void setAbsPos2(int newPos2) {
		absPos2 = newPos2;
	}
	
	public void setProtName2(String newProtName2) {
		protName2=newProtName2;
	}

	public void setProtName1(String newProtName1) {
		protName1=newProtName1;
	}

}
