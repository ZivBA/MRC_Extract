package meshi.applications.TriC.yeast;


import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.applications.TriC.TricYeastAlignment;
import meshi.applications.prediction.GDTcalculator;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Residue;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;

public class PlaceADP  extends MeshiProgram implements Residues,AtomTypes {
	
	public static void main(String[] args) throws Exception {
		init(args);
		
		// Defining the arrangement 
		String topTrue =   "BDAGZQHE"; // OMS
		String chainsTop = "FEAGDHCB";

//		String botTrue =   "BDAGZQHE"; // OMS
//		String chainsBot = "NMIOLPKJ";

		TricYeastAlignment alignments = new TricYeastAlignment();
		AtomList oldList = new AtomList("C:\\Users\\Nir\\TRiC\\Crystallography\\refine_17_pdb_hkl\\refine_17.pdb");
		BufferedWriter bw = new BufferedWriter(new FileWriter("refine_17_fixed_ADP_K.pdb"));
		// Doing top ring first
		for (int chain=0 ; chain<8 ; chain++) {
			AtomList tmpChainList = new AtomList();
			AtomList chainAtoms = oldList.chainFilter(chainsTop.charAt(chain)+"");
			String alignment1Q3R = alignments.getAlignment("I");
			String alignmentChain = alignments.getAlignment(topTrue.charAt(chain)+"");
			int res1Q3RCounter=0;
			int resChainCounter=0;
			for (int res1Q3R=0; res1Q3R<alignment1Q3R.length() ; res1Q3R++) {
				if (alignment1Q3R.charAt(res1Q3R)!='-') {
					res1Q3RCounter++;
				}
				if (alignmentChain.charAt(res1Q3R)!='-') {
					resChainCounter++;
				}
				if ((alignment1Q3R.charAt(res1Q3R)!='-') & (alignmentChain.charAt(res1Q3R)!='-')) {
					Atom atom;
					atom = chainAtoms.findAtomInList("CA", resChainCounter);
					if (atom!=null) {
						Atom newAtom = new Atom(atom.x(), atom.y(), atom.z(),
								"CA", Residue.one2three(alignment1Q3R.charAt(res1Q3R)), res1Q3RCounter, -1);
						newAtom.setChain(""+chainsTop.charAt(chain));
						tmpChainList.add(newAtom);
					}
				}
			}					
			Atom.resetNumberOfAtoms();
			AtomList fragList = new AtomList("1A6E_A_ADP_fragment_withK.pdb");
			GDTcalculator.alignBySubset(tmpChainList, fragList, 0.75);
			fragList.setChain(topTrue.charAt(chain)+"");
			for (int c=0 ; c<fragList.size() ; c++) {
				bw.write(fragList.atomAt(c).toString() + "\n");
			}
			bw.write("TER\n");			
		}
//		// Doing bottom ring second
//		for (int chain=0 ; chain<8 ; chain++) {
//			AtomList tmpChainList = new AtomList();
//			AtomList chainAtoms = oldList.chainFilter(chainsBot.charAt(chain)+"");
//			String alignment1Q3R = alignments.getAlignment("I");
//			String alignmentChain = alignments.getAlignment(botTrue.charAt(chain)+"");
//			int res1Q3RCounter=0;
//			int resChainCounter=0;
//			for (int res1Q3R=0; res1Q3R<alignment1Q3R.length() ; res1Q3R++) {
//				if (alignment1Q3R.charAt(res1Q3R)!='-') {
//					res1Q3RCounter++;
//				}
//				if (alignmentChain.charAt(res1Q3R)!='-') {
//					resChainCounter++;
//				}
//				if ((alignment1Q3R.charAt(res1Q3R)!='-') & (alignmentChain.charAt(res1Q3R)!='-')) {
//					Atom atom;
//					atom = chainAtoms.findAtomInList("CA", resChainCounter);
//					if (atom!=null) {
//						Atom newAtom = new Atom(atom.x(), atom.y(), atom.z(),
//								"CA", Residue.one2three(alignment1Q3R.charAt(res1Q3R)), res1Q3RCounter, -1);
//						newAtom.setChain(""+chainsBot.charAt(chain));
//						tmpChainList.add(newAtom);
//					}
//				}
//			}
//			Atom.resetNumberOfAtoms();
//			AtomList fragList = new AtomList("ADP_frag.pdb");
//			GDTcalculator.alignBySubset(tmpChainList, fragList, 0.75);
//			fragList.setChain(botTrue.toLowerCase().charAt(chain)+"");
//			for (int c=0 ; c<fragList.size() ; c++) {
//				bw.write(fragList.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");			
//		}
		bw.write("END\n");			
		bw.close();		
	}	
	
	protected static void init(String[] args) {
		int zvl = ALA; // force the reading of "meshi.parameters.Residues"
		zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"
		initRandom(333);
	}	

}