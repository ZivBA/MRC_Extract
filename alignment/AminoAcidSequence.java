package alignment;

public class AminoAcidSequence extends Sequence {

	final double gapOpening = -99999999.9;
	final double gapExtension = -0.0; // Note thats this is for a gap against this position.
	
	public AminoAcidSequence(String seq) {
		for (int c=0 ; c<seq.length() ; c++) {
			add(new AminoAcidPosition(seq.charAt(c), gapOpening, gapExtension));
		}
	}

}
