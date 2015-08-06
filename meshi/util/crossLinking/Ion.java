package meshi.util.crossLinking;

public class Ion {

	private double massMonoIsotopic = -1;
	private double intensityMS1 = -1;
	private int charge = -1;
	private int scan = -1;
	private String fileName = null;
	private double retentionTime = -1; // In seconds
	private double[][] ms2Spectrum = null; // size: Nx2 ; first column - M/Z ; second column - intensity
	
	public Ion() {}
	
	public Ion(double massMonoIsotopic, int charge, int scan, String fileName, 
			double intensityMS1, double retentionTime, double[][] ms2Spectrum) {
		this.massMonoIsotopic = massMonoIsotopic;
		this.charge = charge;
		this.scan = scan;
		this.fileName = fileName;
		this.retentionTime = retentionTime;
		this.ms2Spectrum = ms2Spectrum;
		this.setIntensityMS1(intensityMS1);
	}

	public double massMonoIsotopic() {
		return massMonoIsotopic;
	}

	public int charge() {
		return charge;
	}

	public int scan() {
		return scan;
	}

	public String fileName() {
		return fileName;
	}

	public double retentionTime() {
		return retentionTime;
	}

	public double intensityMS1() {
		return intensityMS1;
	}

	public double[][] ms2Spectrum() {
		return ms2Spectrum;
	}

	public void setMassMonoIsotopic(double massMonoIsotopic) {
		this.massMonoIsotopic = massMonoIsotopic;
	}

	public void setCharge(int charge) {
		this.charge = charge;
	}

	public void setScan(int scan) {
		this.scan = scan;
	}

	public void setFileName(String fileName) {
		this.fileName = fileName;
	}

	public void setRetentionTime(double retentionTime) {
		this.retentionTime = retentionTime;
	}

	public void setMs2Spectrum(double[][] ms2Spectrum) {
		this.ms2Spectrum = ms2Spectrum;
	}

	public void setIntensityMS1(double intensityMS1) {
		this.intensityMS1 = intensityMS1;
	}

}
