package meshi.util.crossLinking;

import java.util.Vector;

/**
 * In silico digestion of a single polypeptide
 * @author Nir
 **/
public class Digestion {
	
	private String sequence = null;
	private Peptidase peptidase = null;
	private int[][] peptides = null;
	private int misCleavages = -1;
	private int minLength = -1;

	public Digestion(FastaSeq fasta, Peptidase peptidase, int misCleavages, int minLength) {
		sequence = fasta.seq();
		this.peptidase = peptidase;
		this.misCleavages = misCleavages;
		this.minLength = minLength;
		peptides = digest();
	}
	
	private int[][] digest() {
		Vector<int[]> peps = new Vector<int[]>();
		int startInd = 1;
		int lastInd = 1;
		do {
			lastInd=startInd;
			for (int missedSite = 0 ; (missedSite<=misCleavages) && (lastInd<=sequence.length()) ; missedSite++ ) {
				while (!peptidase.cutC(sequence, lastInd)) {
					lastInd++;
				}
				if (lastInd-startInd+1 >= minLength) {
					int [] tmp = {startInd , lastInd};
					peps.add(tmp);
				}
				lastInd++;				
			}
			while (!peptidase.cutC(sequence, startInd)) {
				startInd++;
			}
			startInd++;
		} while (startInd<=sequence.length());
		
		
		int[][] out = new int[peps.size()][2];
		for (int c=0 ; c<peps.size() ; c++) {
			out[c][0]=peps.get(c)[0];
			out[c][1]=peps.get(c)[1];
 		}
		return out;
	}
	
//	public String toString() {
//		String out = "";
//		for (int c=0 ; c<peptides.length ; c++)
//			out += ((c+1) + ". " + peptides[c][0] + "-" + sequence.substring(peptides[c][0]-1,peptides[c][1]) + 
//					"-" + peptides[c][1] + "\n");
//		return out;
//	}

	public String toString() {
		String out = "";
		for (int c=0 ; c<peptides.length ; c++)
			out += ("'"+ sequence.substring(peptides[c][0]-1,peptides[c][1]) + 
					"',...\n");
		return out;
	}

	public String toString(int pepInd) {
		String out = (peptides[pepInd][0] + "-" + sequence.substring(peptides[pepInd][0]-1,peptides[pepInd][1]) + 
				"-" + peptides[pepInd][1] + "\n");
		return out;
	}
	
	public int[] getLimits(int pepInd) {
		return peptides[pepInd];
	}
	
	public String getSeq(int pepInd) {
		return sequence.substring(peptides[pepInd][0]-1,peptides[pepInd][1]);
	}
	
	public int getNumOfPeptides() {
		return peptides.length;
	}
	
	
	public static void main(String[] str) {
		Digestion digest = new Digestion(new FastaSeq(str[0].trim()), new Trypsin(), 4 , 4);

		System.out.println(digest.getSeq(57) + " " + digest.getSeq(144));
//		System.out.print(digest);

//		MassList massList = new MassList();
//		massList.printMassList(digest);
		
//		MassList massList = new MassList();
//		massList.printXLMassList(digest, digest, str[2].trim());
		
//		MGF mgf = new MGF(str[1].trim());
//		System.out.print(mgf);
	}

}
