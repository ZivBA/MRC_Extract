package meshi.molecularElements.hydrogenBonds.newParameters;


/**
* See what these constants means in heading of the DahiyatHydrogenBond class.
**/
public abstract class DahiyatImplementationConstants {
    // Angle-dependance parameters
	protected abstract double sigmoidBeginsAngleCenterOnH();
	protected abstract double sigmoidEndsAngleCenterOnH();
	protected abstract double sigmoidBeginsNoH();	
	protected abstract double sigmoidEndsNoH();

	// Distance-dependance parameters.
	protected abstract double continueAfterSigmoid();
	protected abstract double valAtp1();
	protected abstract double valAtp2();	
}
