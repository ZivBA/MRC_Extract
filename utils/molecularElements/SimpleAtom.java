package utils.molecularElements;

import static utils.fileUtilities.FileProcessor.*;

/**
 * Created by Ziv_BA on 30/07/2015.
 */
public class SimpleAtom {

	protected String originalString;
	protected short number;
	protected String name;
	protected String aAcidName;
	protected char chain;
	protected short aAcidSequence;
	public SimpleAtom(String atom) {

		originalString = atom;
		name = atom.substring(ATOM_NAME_START, ATOM_NAME_END + 1);
		number = Short.parseShort(atom.substring(ATOM_NUM_START, ATOM_NUM_END + 1).trim());
		aAcidName = atom.substring(RES_NAME_START, RES_NAME_END + 1);
		chain = atom.charAt(CHAIN_ID);
		aAcidSequence = Short.parseShort(atom.substring(RES_SEQ_START, RES_SEQ_END + 1).trim());
	}

	public String getaAcidName() {
		return aAcidName;
	}

	public void setaAcidName(String newAcid){
		this.aAcidName = newAcid;
		originalString = originalString.substring(0,RES_NAME_START)+ newAcid +
				originalString.substring(RES_NAME_END+1);
	}

	public String getOriginalString() {
		return originalString;
	}
}
