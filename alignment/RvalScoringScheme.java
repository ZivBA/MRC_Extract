package alignment;

public class RvalScoringScheme implements ScoringScheme {
	
	public RvalScoringScheme() {	}
	
	public double score(Position pos1, Position pos2) {
		return ((RvalPosition) pos1).fitToAA( ((AminoAcidPosition) pos2).getAminoTypeIndex() );
	}
	
}
