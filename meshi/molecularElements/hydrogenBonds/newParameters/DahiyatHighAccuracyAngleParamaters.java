package meshi.molecularElements.hydrogenBonds.newParameters;

public class DahiyatHighAccuracyAngleParamaters extends
		DahiyatImplementationConstants {

	public final  double sigmoidBeginsAngleCenterOnH = 90.0;
    public final  double sigmoidEndsAngleCenterOnH = 150.0;
    public final  double sigmoidBeginsNoH = 90.0;
    public final  double sigmoidEndsNoH = 110.0;
    public final  double continueAfterSigmoid = 0.6; // Angstroms
    public final  double valAtp1 = 0.98;
    public final  double valAtp2 = 0.04;    

    protected double sigmoidBeginsAngleCenterOnH() {
		return sigmoidBeginsAngleCenterOnH;
	}

	protected double sigmoidBeginsNoH() {
		return sigmoidBeginsNoH;
	}

	protected double sigmoidEndsAngleCenterOnH() {
		return sigmoidEndsAngleCenterOnH;
	}

	protected double sigmoidEndsNoH() {
		return sigmoidEndsNoH;
	}

	protected double continueAfterSigmoid() {
		return continueAfterSigmoid;
	}
	
	protected double valAtp1() {
		return valAtp1;
	}
	
	protected double valAtp2() {
		return valAtp2;
	}	
	
	
}
