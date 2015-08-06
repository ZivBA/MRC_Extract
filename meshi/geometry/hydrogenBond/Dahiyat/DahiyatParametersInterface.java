package meshi.geometry.hydrogenBond.Dahiyat;

/**
* Implementing this interface should give parameters constants to various forms of the Dahiyat formulation.
**/
public interface DahiyatParametersInterface {
	// Donor-Acceptor dependence parameters. Each of these methods should return a matrix that specify
	// some Donor-Acceptor property in the HB between Tsai atom types. The matrices are therefore 14x14.
	// The properties are related to the sigmoid shape of "HB_Sigma".
	public double[][] end();
	public double[][] start();	
	public double[][] p1();	
	public double[][] p2();	
	public double[][] valAtp1();	
	public double[][] valAtp2();
	
	// Angle dependence parameters.
	// ----------------------------
	// When the middle atom of the angle is a hydrogen 
	public double sigmoidBeginsWithH();
	public double sigmoidEndsWithH();
	// When the middle atom of the angle is NOT a hydrogen (meaning no hydrogen in the angle at all) 
	public double sigmoidBeginsNoH();	
	public double sigmoidEndsNoH();
}
