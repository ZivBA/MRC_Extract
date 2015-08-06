package meshi.molecularElements.hydrogenBonds;

public class DahiyatHighAccuracyAngleParamaters extends
		DahiyatImplementationConstants {

	public final  double sigmoidBeginsAngleCenterOnH = 90.0;
    public final  double sigmoidEndsAngleCenterOnH = 150.0;
    public final  double sigmoidBeginsNoH = 90.0;
    public final  double sigmoidEndsNoH = 110.0;
	
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

}
