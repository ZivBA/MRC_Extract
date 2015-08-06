package meshi.applications.Valpuesta;

import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.overlap.Overlap;


public class FitTMbackTo2XSM extends MeshiProgram implements Residues,AtomTypes  {

	public static void main(String[] args) {
		init(args);

		String outfile = "2XSM_as_1Q3R.pdb";
		String chains="ABCDEFGHIJKLMNOP";
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
			for (int let=0 ; let<chains.length() ; let++) {
				Atom.resetNumberOfAtoms();
				AtomList out = alignTMback(""+chains.charAt(let));
				for (int c=0 ; c<out.size() ; c++) {
					bw.write(out.atomAt(c).toString() + "\n");
				}
				bw.write("TER\n");			
			}
			bw.write("END\n");			
			bw.close();
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}
	}


	protected static AtomList alignTMback(String let) {
		AtomList ref = (new AtomList("2XSM_allALA.pdb")).chainFilter(let);
		boolean doEQ=true;
		boolean doMI=true;
		boolean doAP=true;		
		if (let.equals("I") | let.equals("J") | let.equals("L") | let.equals("M") | let.equals("N") | let.equals("O")) {
			doMI = false;
			doAP = false;
		}		
		if (let.equals("K")) {
			ref = (new AtomList("2XSM_allALA_full_K_for_TM_align.pdb")).chainFilter(let);
			doAP = false;
		}		
		int[] statPartListRes = new int[ref.size()];
		for (int c=0 ; c<ref.size() ; c++) {
			statPartListRes[c] = ref.atomAt(c).residueNumber();
		}
		AtomList out = new AtomList();
		AtomList moveEQ = null;
		AtomList moveMI = null;
		AtomList moveAP = null;
		if (doEQ) {
			moveEQ = new AtomList("TMalignment/"+let+"_EQ.pdb");
			System.out.println(Overlap.fitProteins(ref, statPartListRes, moveEQ, statPartListRes));
		}
		if (doMI) {
			moveMI = new AtomList("TMalignment/"+let+"_MI.pdb");
			System.out.println(Overlap.fitProteins(ref, statPartListRes, moveMI, statPartListRes));
		}
		if (doAP) {
			moveAP = new AtomList("TMalignment/"+let+"_AP.pdb");
			System.out.println(Overlap.fitProteins(ref, statPartListRes, moveAP, statPartListRes));
		}
		for (int c=0 ; c<=149 ; c++) {
			Atom atom = moveEQ.findAtomInListReturningAtom("CA", "B", c);
			if (atom != null) {
				out.add(atom);
			}
		}
		if (doMI) {
			for (int c=150 ; c<=217 ; c++) {
				Atom atom = moveMI.findAtomInListReturningAtom("CA", "B", c);
				if (atom != null) {
					out.add(atom);
				}
			}
		}
		if (doAP) {
			for (int c=218 ; c<=368 ; c++) {
				Atom atom = moveAP.findAtomInListReturningAtom("CA", "B", c);
				if (atom != null) {
					out.add(atom);
				}
			}
		}
		if (doMI) {
			for (int c=369 ; c<=403 ; c++) {
				Atom atom = moveMI.findAtomInListReturningAtom("CA", "B", c);
				if (atom != null) {
					out.add(atom);
				}
			}
		}
		for (int c=404 ; c<=600 ; c++) {
			Atom atom = moveEQ.findAtomInListReturningAtom("CA", "B", c);
			if (atom != null) {
				out.add(atom);
			}
		}
		out.setChain(let);
		return out;
	}

	
	protected static void init(String[] args) {
		int zvl = ALA; // force the reading of "meshi.parameters.Residues"
		zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"
		initRandom(333);
	}	

	
	
}
