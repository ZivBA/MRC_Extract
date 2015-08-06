package meshi.util.crossLinking;

import meshi.molecularElements.Protein;

public class MStools implements ResidueMasses {
	
	
	/**
	 * Gives an 2-array [start,end] with the start and end residue numbers of the trypsin cut 
	 * around the modified residue number 'n' in protein 'prot'   
	 */
	public static int[] trypsinCutAroundModifiedResidue(Protein prot, int n) {
		Trypsin trypsin = new Trypsin();
		int start = n;
		while (!trypsin.cutN(prot.getSequence(), start)) {
			start--;
		}
		int end = n+1;
		while (!trypsin.cutC(prot.getSequence(), end)) {
			end++;
		}
		int[] tmp = {start,end};
		return tmp;
	}
	

	/**
	 * Gives the mass of the peptide who's sequence is given in 'seq'.
	 **/
	public static double massOfPeptide(int[] seq) {
		double mass = 0.0;
		for (int c=0 ; c<seq.length ; c++) {
			mass += MW_AA[seq[c]];
		}
		return mass;
	}	
	
}
