package meshi.util.crossLinking;

public class Peptide {
	
	private String sequence = null;
	private int posStartInProt = -1;
	private int posEndInProt = -1;

	public Peptide(String sequence, int posStartInProt, int posEndInProt) {
		this.sequence = sequence;
		this.posStartInProt = posStartInProt;
		this.posEndInProt = posEndInProt;
	}
	
	
	public String sequence() {
		return sequence;
	}

	public int posStartInProt() {
		return posStartInProt;
	}

	public int posEndInProt() {
		return posEndInProt;
	}


}
