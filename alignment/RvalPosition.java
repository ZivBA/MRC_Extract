package alignment;

public class RvalPosition implements Position {

	private double[] Rfits;
	private int resType;
	private int resNum;
	private double gapOpeningScore =  -999999999.9;
	private double gapAligningScore = -0.25;
	
	public RvalPosition(int resNum, int resType, double[] Rfits) {
		this.Rfits = Rfits;
		this.resNum = resNum;
		this.resType = resType;
	}

	@Override
	public double gapOpeningScore() {
		return gapOpeningScore;
	}

	public void setGapOpeningScore(double newScore) {
		gapOpeningScore = newScore;
	}
	
	@Override
	public double gapAligningScore() {
		return gapAligningScore;
	}

	@Override
	public String string() {
		String aa = "ACDEFGHIKLMNPQRSTVWY";
		return aa.substring(resType, resType+1);
	}

	@Override
	public String gapString() {
		return "-";
	}
	
	public double fitToAA(int aaType) {
		return Rfits[aaType];
	}

}
