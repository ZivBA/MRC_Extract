package meshi.applications.TriC.yeast;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.DecimalFormat;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Residue;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;

public class Willison64possibilities extends MeshiProgram implements Residues,AtomTypes  {

	public static void main(String[] args) {
		init(args);
		
		String[] seqs = File2StringArray.f2a("/Users/nirka/TRiC/New_Work_on_TRiC_14_8_2011/Gunnar/Yeast_Thermosome_alignment_1line");
		String seq_13QR = seqs[0].substring(8);
		String seq_1A6DA = seqs[1].substring(8);
		String seq_1A6DB = seqs[2].substring(8);
		String seq_A = seqs[3].substring(8);
		String seq_B = seqs[4].substring(8);
		String seq_G = seqs[5].substring(8);
		String seq_D = seqs[6].substring(8);
		String seq_E = seqs[7].substring(8);
		String seq_H = seqs[8].substring(8);
		String seq_Q = seqs[9].substring(8);
		String seq_Z = seqs[10].substring(8);
		
/*		// Checking consistency of sequences
		System.out.println("Starting");
		String seq = seq_Q.replaceAll("-", "");
		AtomList al = new AtomList("/Users/nirka/TRiC/New_Work_on_TRiC_14_8_2011/Gunnar/3P9D.pdb").chainFilter("H");
		for (int c=0 ; c<al.size() ; c++) {
			if (al.atomAt(c).name().equals("CA")) {
				if (!Residue.one2three(seq.charAt(al.atomAt(c).residueNumber()-1)).equals(al.atomAt(c).residueName())) {
					System.out.println("Mismatch: " + al.atomAt(c).residueNumber() + " " +
							al.atomAt(c).residueName() + " " + Residue.one2three(seq.charAt(al.atomAt(c).residueNumber()-1)));
				}
			}				
		}
		System.out.println("Ending"); */
// Removed 3 residues with numbering >1000 in chain C,K,c,k at position 370. It's a very strange mistake in the PDB deposit!
// G345 was changed to D because Willison used this induced mutation.
		
		
		// Finding the positions to work with:
		boolean[] toTake = new boolean[seq_A.length()];
		for (int c=0 ; c<toTake.length ; c++) {
			toTake[c] = true;
			for (int d=3 ; d<=10 ; d++) { // Any gap invalidates
				if (seqs[d].charAt(c)=='-') {
					toTake[c] = false;
				}
			}
		}
		String chains = "ABCDEFGH";
		int[] chains2genes = {3,4,5,6,7,10,8,9};
		for (int d=0 ; d<chains.length() ; d++) {
			AtomList al = new AtomList("/Users/nirka/TRiC/New_Work_on_TRiC_14_8_2011/Gunnar/3P9D.pdb").chainFilter(""+chains.charAt(d));
			String seq = seqs[chains2genes[d]];
			int nonGapCounter = 0;
			for (int c=0 ; c<toTake.length ; c++) {
				if (seq.charAt(c)!='-') {
					nonGapCounter++;
					if (al.findAtomInList("CA", nonGapCounter)==null) {
						toTake[c] = false;
					}
				}
			}
		}
//		String genes = "ABGDEHQZ";
//		int nonGapCounter = 0;
//		int toTakeCounter = 0;
//		for (int d=3 ; d<=10 ; d++) { 
//			String seq = seqs[d];
//			for (int c=0 ; c<toTake.length ; c++) {
//				if (seq.charAt(c)!='-') {
//					nonGapCounter++;
//					if (toTake[c]) {
//						toTakeCounter++;
//					}
//				}
//			}
//		}
//		System.out.println(toTakeCounter*100.0/nonGapCounter + "% taken.");
//		// On average about 81.7% of the amino acids are represented in the models
						
		
		// Defining the arrangement 
		String topTrue = "ZQHEBDAG";
		String chainsTop = "FEAGDHCB";
		String genesTop =  "ZEAHDQGB";
		String botTrue = "ZQHEBDAG";
		String chainsBot = "NMIOLPKJ";
		String genesBot =  "ZEAHDQGB";

		try {
			int modelNumber = 0;
			for (int delta1=0 ; delta1<8 ; delta1++) { // moving the real 2-fold "delta1" units in the top ring order from willisons 2-fold.
				String outputFile = "model_"+modelNumber+"_"+topTrue+"_"+botTrue+"_"+delta1+".pdb";
				BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
				modelNumber++;
				// Doing top ring on 3P9D
				for (int pos1=0 ; pos1<8 ; pos1++) { // Position on Willison's PDB
					char chainID = chainsTop.charAt(pos1);
					char geneInWillison = genesTop.charAt(pos1);
					char toPut = topTrue.charAt((pos1+8-delta1) % 8);
					String querySeq = seqs[getSeqNum(toPut)];
					String templateSeq = seqs[getSeqNum(geneInWillison)];
					AtomList newList = HMunit(querySeq,templateSeq,toTake,""+chainID,"D");
					newList.setChain(""+chainID);
					for (int c=0 ; c<newList.size() ; c++) {
						bw.write(newList.atomAt(c).toString() + "\n");
					}
					bw.write("TER\n");			
				}
				// Doing bottom ring on 3P9D
				for (int pos1=0 ; pos1<8 ; pos1++) { // Position on Willison's PDB
					char chainID = chainsBot.charAt(pos1);
					char geneInWillison = genesBot.charAt(pos1);
					char toPut = botTrue.charAt((pos1+delta1) % 8);
					String querySeq = seqs[getSeqNum(toPut)];
					String templateSeq = seqs[getSeqNum(geneInWillison)];
					AtomList newList = HMunit(querySeq,templateSeq,toTake,""+chainID,"D");
					newList.setChain(""+chainID);
					for (int c=0 ; c<newList.size() ; c++) {
						bw.write(newList.atomAt(c).toString() + "\n");
					}
					bw.write("TER\n");			
				}
				bw.write("END\n");			
				bw.close();
			}		
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}

		
		
	}
	
	
	
	
	protected static AtomList HMunit(String querySeq,String templateSeq,boolean[] toTake,String chainID, String DorE) {
		AtomList al = new AtomList("/Users/nirka/TRiC/New_Work_on_TRiC_14_8_2011/Gunnar/3P9"+DorE+".pdb").chainFilter(chainID);
		Atom.resetNumberOfAtoms();
		AtomList newList = new AtomList();
		int templateCounter = 0;
		int queryCounter = 0;
		for (int resC=0 ; resC<toTake.length ; resC++) {
			if (templateSeq.charAt(resC)!='-') {
				templateCounter++;
			}
			if (querySeq.charAt(resC)!='-') {
				queryCounter++;
			}
			if (toTake[resC]) {
				Atom atom;
				// Doing N
				atom = al.findAtomInList("N", templateCounter);
				newList.add(new Atom(atom.x(), atom.y(), atom.z(),
						   "N", Residue.one2three(querySeq.charAt(resC)), queryCounter, -1));
				// Doing CA
				atom = al.findAtomInList("CA", templateCounter);
				newList.add(new Atom(atom.x(), atom.y(), atom.z(),
						   "CA", Residue.one2three(querySeq.charAt(resC)), queryCounter, -1));
				// Doing C
				atom = al.findAtomInList("C", templateCounter);
				newList.add(new Atom(atom.x(), atom.y(), atom.z(),
						   "C", Residue.one2three(querySeq.charAt(resC)), queryCounter, -1));
				// Doing O
				atom = al.findAtomInList("O", templateCounter);
				newList.add(new Atom(atom.x(), atom.y(), atom.z(),
						   "O", Residue.one2three(querySeq.charAt(resC)), queryCounter, -1));
			}
		}
		return newList;		
	}
	
	protected static int getSeqNum(char gene) {
		switch(gene) {
		case 'A':
			return 3;
		case 'B':
			return 4;
		case 'G':
			return 5;
		case 'D':
			return 6;
		case 'E':
			return 7;
		case 'H':
			return 8;
		case 'Q':
			return 9;
		case 'Z':
			return 10;
		}
		return -1;
	}
	
	
	

	
	
	
	
	protected static void init(String[] args) {
		int zvl = ALA; // force the reading of "meshi.parameters.Residues"
		zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"
		initRandom(333);
	}	

	
	
}
