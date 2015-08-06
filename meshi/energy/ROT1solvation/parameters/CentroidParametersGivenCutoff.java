package meshi.energy.ROT1solvation.parameters;

public class CentroidParametersGivenCutoff extends AbstractCBParameters {

	private double cutoff = -999;
	private double sigmoidRange = 0.1;	
	
	public CentroidParametersGivenCutoff(String path, String cutoffString) {
		super(path + "/CENTROID_" + cutoffString + ".txt");
		cutoff = (new Double(cutoffString)).doubleValue();
	}
	
	protected double cutoff() {
		return cutoff;
	}
	
	protected double sigmoidRange() {
		return sigmoidRange;
	}
	
}