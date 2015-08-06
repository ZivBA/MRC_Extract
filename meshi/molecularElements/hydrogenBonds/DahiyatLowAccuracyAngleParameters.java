package meshi.molecularElements.hydrogenBonds;

/**
 * 
 * This class is mainly good for side-chain modeling, where the rotamer approximation
 * compromises the angular accuracy of the HB.
 * 
 **/
public class DahiyatLowAccuracyAngleParameters extends
		DahiyatImplementationConstants {

	public final  double sigmoidBeginsAngleCenterOnH = 90.0;
    public final  double sigmoidEndsAngleCenterOnH = 110.0;
    public final  double sigmoidBeginsNoH = 80.0;
    public final  double sigmoidEndsNoH = 100.0;
	
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
