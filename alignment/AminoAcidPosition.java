package alignment;

public class AminoAcidPosition implements Position {

	protected String aminoAcidLetter = "";
	protected int aminoAcidInd = -1;
	protected double gapOpeningScore;
	protected double gapAligningScore;
	
	public AminoAcidPosition(char aminoAcid , double gapOpeningScore , double gapAligningScore) {
		aminoAcidLetter = ""+aminoAcid;
		switch (aminoAcid) {
		case 'A':	aminoAcidInd = 0;
		break;
		case 'C':	aminoAcidInd = 1;
		break;
		case 'D':	aminoAcidInd = 2;
		break;
		case 'E':	aminoAcidInd = 3;
		break;
		case 'F':	aminoAcidInd = 4;
		break;
		case 'G':	aminoAcidInd = 5;
		break;
		case 'H':	aminoAcidInd = 6;
		break;
		case 'I':	aminoAcidInd = 7;
		break;
		case 'K':	aminoAcidInd = 8;
		break;
		case 'L':	aminoAcidInd = 9;
		break;
		case 'M':	aminoAcidInd = 10;
		break;
		case 'N':	aminoAcidInd = 11;
		break;
		case 'P':	aminoAcidInd = 12;
		break;
		case 'Q':	aminoAcidInd = 13;
		break;
		case 'R':	aminoAcidInd = 14;
		break;
		case 'S':	aminoAcidInd = 15;
		break;
		case 'T':	aminoAcidInd = 16;
		break;
		case 'V':	aminoAcidInd = 17;
		break;
		case 'W':	aminoAcidInd = 18;
		break;
		case 'Y':	aminoAcidInd = 19;
		break;

		// special cases
		case 'U':   aminoAcidInd = 1;
		break;


		default: throw new RuntimeException("Unkown amino acid type: " + aminoAcid);
		}
		this.gapOpeningScore = gapOpeningScore;
		this.gapAligningScore = gapAligningScore;
	}

	@Override
	public double gapOpeningScore() {
		return gapOpeningScore;
	}

	@Override
	public double gapAligningScore() {
		// TODO Auto-generated method stub
		return gapAligningScore;
	}

	@Override
	public String string() {
		return aminoAcidLetter;
	}

	@Override
	public String gapString() {
		return "-";
	}
	
	public int getAminoTypeIndex() {
		return aminoAcidInd;
	}

}
