package meshi.energy.ROT1solvation.parameters;

public class CBParametersGivenCutoff extends AbstractCBParameters {

	private double cutoff = -999;
	private double sigmoidRange = 0.1;	
	
	public CBParametersGivenCutoff(String path, String cutoffString) {
		super(path + "/CB_" + cutoffString + ".txt");
		cutoff = (new Double(cutoffString)).doubleValue();
	}
	
	protected double cutoff() {
		return cutoff;
	}
	
	protected double sigmoidRange() {
		return sigmoidRange;
	}
	
}
