package meshi.util.external;

import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.applications.TriC.TricAlignment;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Residue;
import meshi.parameters.Residues;
import meshi.util.file.MeshiWriter;

public class PrepareSCWRL implements Residues {

	private AtomList al = null;
	
	public PrepareSCWRL(AtomList al) {
		this.al = al;		
	}
	
	public String getSCWRLsequence() {
		int zvl = ALA;
		String output = "";
		for (int c=0 ; c<al.size() ; c++) {
			if (al.atomAt(c).name().equals("CA")) {
				output += Residue.three2one(al.atomAt(c).residueName());
			}
		}
		return output;		
	}

/**
 * This method will return a SCWRL sequence file, which will free for side-chain
 * modeling only residues with Ca's that are less than 'disFromAnotherChain' from 
 * another chain.
 **/
	public String getSCWRLsequenceForInterface(double disFromAnotherChain) {
		int zvl = ALA;
		String output = "";
		al.freeze();
		for (int c1=0 ; c1<al.size() ; c1++) {
			for (int c2=c1+1 ; c2<al.size() ; c2++) {
				if (al.atomAt(c1).name().equals("CA") && (al.atomAt(c2).name().equals("CA")) &&
						!al.atomAt(c2).chain().equals(al.atomAt(c1).chain())) {
					if (al.atomAt(c1).distanceFrom(al.atomAt(c2)) < disFromAnotherChain) {
						al.atomAt(c1).defrost();
						al.atomAt(c2).defrost();
					}							
				}				
			}
		}
		boolean lastIsChainA = true;
		for (int c=0 ; c<al.size() ; c++) {
			if (al.atomAt(c).name().equals("CA")) {
				if (al.atomAt(c).chain().equals("B")) {
					if (lastIsChainA) {
						output += "\n";
					}
					lastIsChainA = false;
				}
				if (al.atomAt(c).frozen()) {
					output += Residue.three2one(al.atomAt(c).residueName()).toLowerCase();
				}
				else {
					output += Residue.three2one(al.atomAt(c).residueName());
				}
			}
		}
		return output;		
	}
	
	/**
	 * This method will return a SCWRL sequence file, which will free for side-chain
	 * modeling only residues that are not conserved in the alignment.
	 * Assuming template sequence starts at residue 1.
	 **/
		public String getSCWRLsequenceForAlignment(String queryAlignment,String templateAlignment) {
			String output = "";
			int resCounter = 1;
			int iden = 0;
			for (int c=0 ; c<queryAlignment.length() ; c++) {
				if ((queryAlignment.charAt(c)!='-') && (templateAlignment.charAt(c)!='-') &&
						(al.findAtomInList("CA", resCounter)!=null)) {
					if (queryAlignment.charAt(c)==templateAlignment.charAt(c)) {
						output += ("" + queryAlignment.charAt(c)).toLowerCase();
						iden++;
					}
					else {
						output += queryAlignment.charAt(c);
					}
				}
				if (templateAlignment.charAt(c)!='-') {
					resCounter++;
				}
			}
			System.out.println("Sequence identity: " + iden*100.0/output.length() + "%");
			return output;		
		}

	
		public static void main(String[] args) {
			TricAlignment tric = new TricAlignment();
			String templateS = tric.getAlignment("I");
			String queryS = tric.getAlignment(args[0]);
			PrepareSCWRL prepSCWRL = new PrepareSCWRL(new AtomList(args[1]));
			String scwrlString = prepSCWRL.getSCWRLsequenceForAlignment(queryS, templateS);
			try{
				BufferedWriter bw = new BufferedWriter(new FileWriter(args[2]));
				for (int c=0 ; c<scwrlString.length() ; c++) {
					if (c!=0 && ((c/50==c/50.0))) {
						bw.write("\n");
						System.out.println();
					}
					System.out.print(scwrlString.charAt(c));
					bw.write(scwrlString.charAt(c));
				}
				bw.close();
			}
			catch(Exception e) {
				throw new RuntimeException(e.getMessage());
			}    			
			System.out.println();
		}
		
//	public static void main(String[] args) {
//		PrepareSCWRL prepSCWRL = new PrepareSCWRL(new AtomList(args[0]));
//		String scwrlString = prepSCWRL.getSCWRLsequenceForInterface(12.0);
//		try{
//			BufferedWriter bw = new BufferedWriter(new FileWriter(args[1]));
//			for (int c=0 ; c<scwrlString.length() ; c++) {
//				if (c!=0 && ((c/50==c/50.0))) {
//					bw.write("\n");
//					System.out.println();
//				}
//				System.out.print(scwrlString.charAt(c));
//				bw.write(scwrlString.charAt(c));
//			}
//		    bw.close();
//		}
//		catch(Exception e) {
//		    throw new RuntimeException(e.getMessage());
//		}    			
//		System.out.println();
//	}
}

