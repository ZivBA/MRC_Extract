package meshi.util.crossLinking;

import java.text.DecimalFormat;

import meshi.molecularElements.Atom;
import meshi.molecularElements.Protein;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.KeyWords;
import meshi.util.dssp.DSSP;


public class crossLinkingTools implements Residues, AtomTypes , KeyWords, ResidueMasses {
	
	public static void findCandidates(Protein prot1, Protein prot2, DSSP dssp, String chain1, String chain2, double lowDisTH, double highDisTH, double accessibilityTH) { 
		int NallCL = 0;
		int NaccessCL = 0;
		DecimalFormat fmt = new DecimalFormat("0.##");
		DecimalFormat fmt4 = new DecimalFormat("0.####");
		System.out.println(chain1 + " + " + chain2 + "\n------\n");
		for (int c1=0 ; c1<prot1.atoms().size() ; c1++) {
			for (int c2=0 ; c2<prot2.atoms().size() ; c2++) {
				Atom a1 = prot1.atoms().atomAt(c1);
				Atom a2 = prot2.atoms().atomAt(c2);
				if ((a1.chain().equals(chain1) && a2.chain().equals(chain2)) || 
						(a1.chain().equals(chain2) && a2.chain().equals(chain1))) {
					if ((((a1.type == atomLinked1) && (a2.type == atomLinked2)) || ((a1.type == atomLinked2) && (a2.type == atomLinked1)))) {
						double distance = a1.distanceFrom(a2);
						if ((distance>lowDisTH) & (distance<highDisTH)) {
							NallCL++;
							if ((dssp.relACCofRes(a1.residueNumber(), a1.chain().charAt(0))>accessibilityTH) &&
									(dssp.relACCofRes(a2.residueNumber(), a2.chain().charAt(0))>accessibilityTH)) {
								NaccessCL++;
								System.out.println("CHAIN: " + a1.chain() + " K" + a1.residueNumber() + " <-- " + 
										fmt.format(distance) + " --> " + "CHAIN: " + a2.chain() + " K" + a2.residueNumber() + 
										"      Solv Acc:" + dssp.relACCofRes(a1.residueNumber(), a1.chain().charAt(0)) +
										"  " + dssp.relACCofRes(a2.residueNumber(), a2.chain().charAt(0)) + 
										"    Mass: " + fmt4.format(massOfCrossLinkedPeptide(prot1, a1.residueNumber(), prot2, a2.residueNumber())) +
										"    Score: " + fmt.format(scoreCrossLink(prot1, a1.residueNumber(), prot2, a2.residueNumber(), dssp)));
							}
						}
					}
				}
			}
		}
		System.out.println("------------------------------------\nAll possible: " + NallCL);
		System.out.println("Good SolvACC: " + NaccessCL + " (" + fmt.format(NaccessCL*100.0/NallCL) + "%)");

	}

	public static double massOfCrossLinkedPeptide(Protein prot1, int resCLinProt1, Protein prot2, int resCLinProt2) {
		int[] firstPeptide = MStools.trypsinCutAroundModifiedResidue(prot1, resCLinProt1);
		int[] secondPeptide = MStools.trypsinCutAroundModifiedResidue(prot2, resCLinProt2);
		int[] seq1 = Protein.getSeqOfProt(prot1, firstPeptide[0], firstPeptide[1]);
		int[] seq2 = Protein.getSeqOfProt(prot2, secondPeptide[0], secondPeptide[1]);
		double totalMass = MW_CL + MStools.massOfPeptide(seq1) + MStools.massOfPeptide(seq2) + 2*MW_OH + 2*MW_H;
		return totalMass;
	}
	
	
	public static double scoreCrossLink(Protein prot1, int resCLinProt1, Protein prot2, int resCLinProt2, DSSP dssp) {
		double score = 0;
		Trypsin tripsin = new Trypsin();
		int[] firstPeptide = MStools.trypsinCutAroundModifiedResidue(prot1, resCLinProt1);
		int[] secondPeptide = MStools.trypsinCutAroundModifiedResidue(prot2, resCLinProt2);
		int[] seq1 = Protein.getSeqOfProt(prot1, firstPeptide[0], firstPeptide[1]);
		int[] seq2 = Protein.getSeqOfProt(prot2, secondPeptide[0], secondPeptide[1]);
		if ((seq1.length<3) || (seq2.length<3)) // Removing KK or KR links
			return 66;
		if (!tripsin.isCterm(prot1, firstPeptide[1]))
			if ((prot1.residue(firstPeptide[1]).type==8) && (dssp.relACCofRes(firstPeptide[1],prot1.residue(firstPeptide[1]).ca().chain().charAt(0))>0.15))
				score++;
		if (!tripsin.isCterm(prot1, firstPeptide[1]))
			if ((prot2.residue(secondPeptide[1]).type==8) && (dssp.relACCofRes(secondPeptide[1],prot2.residue(secondPeptide[1]).ca().chain().charAt(0))>0.15))
				score++;
		if (!tripsin.isNterm(prot1, firstPeptide[0]))
			if ((prot1.residue(firstPeptide[0]-1).type==8) && (dssp.relACCofRes(firstPeptide[0]-1,prot1.residue(firstPeptide[0]).ca().chain().charAt(0))>0.15))
				score++;
		if (!tripsin.isNterm(prot2, secondPeptide[0]))
			if ((prot2.residue(secondPeptide[0]-1).type==8) && (dssp.relACCofRes(secondPeptide[0]-1,prot2.residue(secondPeptide[0]).ca().chain().charAt(0))>0.15))
				score++;
		return score;
	}


}
