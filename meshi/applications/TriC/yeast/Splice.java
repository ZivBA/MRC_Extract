package meshi.applications.TriC.yeast;

import java.io.IOException;

import meshi.applications.prediction.GDTcalculator;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.file.MeshiWriter;

public class Splice {

	/**
	 * The purpose of this class is to put changes made into one subunit in one ring into all the 4 occurrences of the subunit in the unit cell.
	 */
	
	
	public static AtomList splice(AtomList oldList , AtomList newList , int firstResToTake , int lastResToTake, String chainName) {
		System.out.println("\n\n" + oldList.atomAt(0) + "\n" + newList.atomAt(0) + "\n" + chainName);
		AtomList outList = new AtomList();
		AtomList oldListUnchanged = new AtomList();
		for (int c=0 ; c<oldList.size() ; c++) {
			if ((oldList.atomAt(c).residueNumber()<firstResToTake) | (oldList.atomAt(c).residueNumber()>lastResToTake)) {
				oldListUnchanged.add(oldList.atomAt(c));
			}
		}
		GDTcalculator.alignBySubset(oldListUnchanged, newList, 0.05);
		boolean put = false;
		for (int c=0 ; c<oldList.size() ; c++) {
			if ((oldList.atomAt(c).residueNumber()>=firstResToTake) & (oldList.atomAt(c).residueNumber()<=lastResToTake)) {
				if (!put) {
					put = true;
					for (int d=0 ; d<newList.size(); d++) {
						if ((newList.atomAt(d).residueNumber()>=firstResToTake) & (newList.atomAt(d).residueNumber()<=lastResToTake)) {
							Atom atom = new Atom(newList.atomAt(d));
							atom.setChain(chainName);
							outList.add(atom);
						}
					}
				}
			}
			else {
				outList.add(oldList.atomAt(c));
			}					
		}
		return outList;
	}
	

	public static AtomList fixAtomNumbers(AtomList oldList) {
		AtomList outList = new AtomList();

		String oldChain = "-";
		for (int c=0 ; c<oldList.size() ; c++) {
			if (oldList.atomAt(c).chain().equals("F") & oldChain.equals("-")) {
				Atom.resetNumberOfAtoms();
				System.out.println(new Atom(0, 0, 0, "X", "X", 0, 0));
				System.out.println("Reseted the atom Number");
			}
			if (oldList.atomAt(c).chain().equals("N") & oldChain.equals("B")) {
				Atom.resetNumberOfAtoms();
				System.out.println("Reseted the atom Number");
				System.out.println(new Atom(0, 0, 0, "X", "X", 0, 0));
			}
			if (oldList.atomAt(c).chain().equals("f") & oldChain.equals("J")) {
				Atom.resetNumberOfAtoms();
				System.out.println("Reseted the atom Number");
				System.out.println(new Atom(0, 0, 0, "X", "X", 0, 0));
			}
			if (oldList.atomAt(c).chain().equals("n") & oldChain.equals("b")) {
				Atom.resetNumberOfAtoms();
				System.out.println("Reseted the atom Number");
				System.out.println(new Atom(0, 0, 0, "X", "X", 0, 0));
			}
			Atom atom = new Atom(oldList.atomAt(c));
			if (atom.name().charAt(0)!='H') {
				outList.add(atom);
			}
			oldChain = atom.chain().trim();
		}
		return outList;
	}

	
	public static void main(String[] args) {
		String fileOld = args[0].trim();
		String fileNew = args[1].trim();
		String fileWrite = args[2].trim();
		String affectedChain = args[3].trim();
		int    firstRes = Integer.parseInt(args[4].trim());
		int    lastRes = Integer.parseInt(args[5].trim());
		
		AtomList oldList = new AtomList(fileOld);
		AtomList newList = new AtomList(fileNew);
		String chainTop1 = "FEAGDHCB";
		String chainBot1 = "NMIOLPKJ";
		String chainTop2 = "feagdhcb";
		String chainBot2 = "nmiolpkj";
		
		int indOfChain = chainTop1.indexOf(affectedChain);
		
		AtomList writeList = new AtomList();
		for (int c=0 ; c<chainTop1.length() ; c++) {
			if (c==indOfChain) {
				AtomList changedList = splice(oldList.chainFilter(""+chainTop1.charAt(c)) , newList.chainFilter(affectedChain) , firstRes , lastRes, ""+chainTop1.charAt(c));
				writeList.add(changedList);
			}
			else {
				writeList.add(oldList.chainFilter(""+chainTop1.charAt(c)));
			}
		}
		for (int c=0 ; c<chainBot1.length() ; c++) {
			if (c==indOfChain) {
				AtomList changedList = splice(oldList.chainFilter(""+chainBot1.charAt(c)) , newList.chainFilter(affectedChain) , firstRes , lastRes, ""+chainBot1.charAt(c));
				writeList.add(changedList);
			}
			else {
				writeList.add(oldList.chainFilter(""+chainBot1.charAt(c)));
			}
		}
		for (int c=0 ; c<chainTop2.length() ; c++) {
			if (c==indOfChain) {
				AtomList changedList = splice(oldList.chainFilter(""+chainTop2.charAt(c)) , newList.chainFilter(affectedChain) , firstRes , lastRes, ""+chainTop2.charAt(c));
				writeList.add(changedList);
			}
			else {
				writeList.add(oldList.chainFilter(""+chainTop2.charAt(c)));
			}
		}
		for (int c=0 ; c<chainBot2.length() ; c++) {
			if (c==indOfChain) {
				AtomList changedList = splice(oldList.chainFilter(""+chainBot2.charAt(c)) , newList.chainFilter(affectedChain) , firstRes , lastRes, ""+chainBot2.charAt(c));
				writeList.add(changedList);
			}
			else {
				writeList.add(oldList.chainFilter(""+chainBot2.charAt(c)));
			}
		}		
		
		AtomList writeListFinal = fixAtomNumbers(writeList);
		try {
			writeListFinal.print(new MeshiWriter(fileWrite));
		} catch (IOException e) {
			e.printStackTrace();
		}
		



	}

}
