package meshi.applications.TriC.yeast;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;

import meshi.applications.TriC.TricYeastAlignment;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;

public class FindContactsForFigure {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		AtomList list = new AtomList("equaZ_equaQ.pdb");
		TricYeastAlignment alignment = new TricYeastAlignment();
		boolean[] interacting = interfacingResidues(list, "Q");
		for (int c=0 ; c<interacting.length ; c++) {
			if (interacting[c]) {
				System.out.print(alignment.getNewResNum('Q', c, 'Q')+",");
			}
		}

	}

	
	
	
	/**
	 * The output is a boolean[2000] array (I assume that the desired chains is not more 1999 residues). True is interfacing.
	 */
	public static boolean[] interfacingResidues(AtomList list, String desiredChain) {
		double CONTACT_TH = 3.8;
		boolean[] out = new boolean[2000];
		for (int c=0 ; c<2000 ; c++) {
			out[c] = false;
		}
		// Looking from A to B
		for (int c=0 ; c<list.size() ; c++) {
			for (int d=0 ; d<list.size() ; d++) {
				if (list.atomAt(c).chain().equals(desiredChain) &&  !list.atomAt(d).chain().equals(desiredChain) &&
						(list.atomAt(c).distanceFrom(list.atomAt(d)) < CONTACT_TH)) {
					out[list.atomAt(c).residueNumber()] = true;
				}
			}
		}
		return out;
	}

}
