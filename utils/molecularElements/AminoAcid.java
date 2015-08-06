package utils.molecularElements;

import java.util.Iterator;
import java.util.List;

/**
 * Created by Ziv_BA on 30/07/2015.
 */
public class AminoAcid implements Iterable<SimpleAtom>{
	private SimpleAtom[] atoms;
	private String name;
	private char chainID;
	private int seqNum;

	/**
	 * constructor for AminoAcid object, gets a string list (such as from PDB file).
	 *
	 * @param listOfAtoms
	 */
	public AminoAcid(List<String> listOfAtoms) {
		atoms = new SimpleAtom[listOfAtoms.size()];

		for (int i = 0; i < listOfAtoms.size(); i++) {
			atoms[i] = new SimpleAtom(listOfAtoms.get(i));
		}
		name = atoms[0].getaAcidName();
		chainID = atoms[0].chain;
		seqNum = atoms[0].aAcidSequence;

	}

	public int getSeqNum() {
		return seqNum;
	}

	public char getChainID() {
		return atoms[0].chain;
	}


	public String getName() {

		return name;
	}

	@Override
	public Iterator<SimpleAtom> iterator() {
		return new Iterator<SimpleAtom>() {
			int counter;
			@Override
			public boolean hasNext() {
				return counter <atoms.length;
			}

			@Override
			public SimpleAtom next() {
				return atoms[counter++];
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}

	public void substituteWith(String newAcid) {
		name = newAcid;
		for (SimpleAtom atom : atoms) {
			atom.setaAcidName(newAcid);
		}
	}
}
