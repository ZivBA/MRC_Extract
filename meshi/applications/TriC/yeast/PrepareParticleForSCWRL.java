package meshi.applications.TriC.yeast;

import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;

public class PrepareParticleForSCWRL {

	public static void main(String[] args) {
		AtomList allAtoms = new AtomList("refine_16-half.ENCAD.pdb");
		Atom.resetNumberOfAtoms();
		writeChain(allAtoms , "A");
		Atom.resetNumberOfAtoms();
		writeChain(allAtoms , "B");
		Atom.resetNumberOfAtoms();
		writeChain(allAtoms , "G");
		Atom.resetNumberOfAtoms();
		writeChain(allAtoms , "D");
		Atom.resetNumberOfAtoms();
		writeChain(allAtoms , "E");
		Atom.resetNumberOfAtoms();
		writeChain(allAtoms , "H");
		Atom.resetNumberOfAtoms();
		writeChain(allAtoms , "Q");
		Atom.resetNumberOfAtoms();
		writeChain(allAtoms , "Z");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "A");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "B");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "C");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "D");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "E");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "F");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "G");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "H");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "I");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "J");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "K");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "L");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "M");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "N");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "O");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "P");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "a");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "b");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "c");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "d");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "e");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "f");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "g");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "h");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "i");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "j");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "k");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "l");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "m");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "n");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "o");
//		Atom.resetNumberOfAtoms();
//		writeChain(allAtoms , "p");
	}

	public static void writeChain(AtomList allAtoms, String chain) {
		
		double takeR = 14.0;
		
		boolean[] take = new boolean[allAtoms.size()];
		for (int c=0 ; c<take.length ; c++) {
			take[c] = false;
		}
		
		AtomList chainAtom = allAtoms.chainFilter(chain);
		for (int c=0 ; c<allAtoms.size() ; c++) {
			for (int d=0 ; d<chainAtom.size() ; d++) {
				if (allAtoms.atomAt(c).distanceFrom(chainAtom.atomAt(d)) < takeR) {
					take[c] = true;
				}
			}
		}
		
		for (int c=0 ; c<allAtoms.size() ; c++) {
			if (!take[c]) {
				for (int d=Math.max(0, c-20) ; d<Math.min(take.length, c+20) ; d++) {
					if (take[d] && (allAtoms.atomAt(d).residueNumber() == allAtoms.atomAt(c).residueNumber()) &&
							(allAtoms.atomAt(d).chain().equals(allAtoms.atomAt(c).chain()))) {
						take[c] = true;
						break;
					}
				}				
			}
		}
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter("do_"+chain+".pdb"));
			for (int c=0 ; c<allAtoms.size() ; c++) {
				if (take[c]) {
					Atom atom = new Atom(allAtoms.atomAt(c));
					if (!atom.residueName().equals("ADB")) {
						bw.write(atom.toString() + "\n");
					}
				}
			}
			bw.close();
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}
	}
	
}
