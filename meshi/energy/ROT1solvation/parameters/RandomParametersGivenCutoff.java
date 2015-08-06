package meshi.energy.ROT1solvation.parameters;

public class RandomParametersGivenCutoff extends AbstractROT1Parameters {

	private double cutoff = -999;
	private double sigmoidRange = 0.1;	
	
	public RandomParametersGivenCutoff(String path, String cutoffString) {
		super(path + "/RANDOM_" + cutoffString + ".txt");
		cutoff = (new Double(cutoffString)).doubleValue();
	}
	
	protected double cutoff() {
		return cutoff;
	}
	
	protected double sigmoidRange() {
		return sigmoidRange;
	}
	
}
