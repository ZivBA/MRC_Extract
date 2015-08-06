package meshi.util.crossLinking;

import meshi.molecularElements.Protein;

public abstract class Peptidase {

	public Peptidase() {}

	boolean isCterm(Protein prot, int residueNumber) {
		if (prot.residue(residueNumber+1) == null) // End of the chain
			return true;
		if (prot.residue(residueNumber+1).ca() == null) // End of the chain
			return true;
		return false;
	}
	
	boolean isNterm(Protein prot, int residueNumber) {
		if (prot.residue(residueNumber-1) == null) // End of the chain
			return true;
		if (prot.residue(residueNumber-1).ca() == null) // End of the chain
			return true;
		return false;
	}	
	
	
	abstract boolean cutN(String seq, int residueNumber);
	abstract boolean cutC(String seq, int residueNumber);
}
