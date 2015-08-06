package meshi.energy.oldSolvate;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.util.mathTools.Sigma;

/**
 * A class to calculate the HB angular score and its derivatives. The two angles in the HB made
 * of 4 atoms are defined as follows:
 *
 *   a1~~~~~~a2
 *             .
 *              .
 *               a3~~~~~~a4
 *
 * ...  -  the HB
 * ~~~  -  covelant bonding
 * Angle 1 - between atoms 1,2,3.
 * Angle 2 - between atoms 2,3,4.
 * There could be two cases:
 * 1) The hydrogen atom in the HB is not known explicitly (NoH). In this case a2 and a3 are the polar atoms, 
 * and a1 and a4 are the heavy atoms to which they are attached (base atoms).
 * 2) The hydrogen atom in the HB is known explicitly (WithH). Either a2 or a3 must be a hydrogen. The other atom 
 * COVELANTLY bonded to the hydrogen must be the polar atom. 
 * 
 * The angular score of each angle is a sigmoid that equals 0.0 bellow a certain threshold (sigmoidBegins) and 1.0 
 * above a certain threshold (sigmoidEnds). Between the thresholds it raises smoothly.
 * The HB angular score is the product of the two angular score of angle 1 and 2.
 *
 * How to use this class:
 * ----------------------
 * 1) The constructor gets the sigmoid thresholds as parameters. The threshold values could be different for the two
 * cases (WithH or NoH). Usually the WithH is more restrictive, and sigmoidBeginsWithH is higher than 
 * sigmoidBeginsNoH.
 * 2) Call the updateAndEvaluateAtoms1234 method giving it the 4 atoms. 
 * 3) The HB angular score and its derivatives with respect to the XYZ coordinates of the four atoms, can now be 
 * accessed through the appropriate getter classes. 
 *
 *
 * The 'zeroDerivative' field:
 * ---------------------------
 * In many cases  (mainly outside the sigmoid transition region [sigmoidBegins,sigmoidEnds]) the derivative
 * of the HB angle score is zero with respect to all the atom coordinates. The  'zeroDerivative' boolean field 
 * is updated to 'true' in those cases, by the in 'updateAndEvaluateAtoms1234' method. If not all the derivatives
 * are zero it is set to 'false'. Using this field is capable of saving a lot of unrelevent operations in 
 * the SolvateEnergy class. 
 **/


public class SolvateHBAngle {
			
	private DistanceMatrix dm;
	private final double sigmoidBeginsWithH;
	private final double sigmoidEndsWithH;
	private final double sigmoidBeginsNoH;
	private final double sigmoidEndsNoH;
	private double sigmoidBegins;
	private double sigmoidEnds;
	private double	hbAngScore;    // This is a variable solvate uses, and has a getter method 
	private boolean	zeroDerivative;    // This is a variable solvate uses, and has a getter method 
	private double	DhbAngScoreDx1; // This is a variable solvate uses, and has a getter method
	private double	DhbAngScoreDy1; // This is a variable solvate uses, and has a getter method
	private double	DhbAngScoreDz1; // This is a variable solvate uses, and has a getter method
	private double	DhbAngScoreDx2; // This is a variable solvate uses, and has a getter method
	private double	DhbAngScoreDy2; // This is a variable solvate uses, and has a getter method
	private double	DhbAngScoreDz2; // This is a variable solvate uses, and has a getter method
	private double	DhbAngScoreDx3; // This is a variable solvate uses, and has a getter method
	private double	DhbAngScoreDy3; // This is a variable solvate uses, and has a getter method
	private double	DhbAngScoreDz3; // This is a variable solvate uses, and has a getter method
	private double	DhbAngScoreDx4; // This is a variable solvate uses, and has a getter method
	private double	DhbAngScoreDy4; // This is a variable solvate uses, and has a getter method
	private double	DhbAngScoreDz4; // This is a variable solvate uses, and has a getter method
	private double	sigmCosAng1;
	private boolean derivativeAng1Zero;
	private double	DsigmCosAng1Dx1;
	private double	DsigmCosAng1Dy1;
	private double	DsigmCosAng1Dz1;
	private double	DsigmCosAng1Dx2;
	private double	DsigmCosAng1Dy2;
	private double	DsigmCosAng1Dz2;
	private double	DsigmCosAng1Dx3;
	private double	DsigmCosAng1Dy3;
	private double	DsigmCosAng1Dz3;
	private double	sigmCosAng2;
	private boolean derivativeAng2Zero;
	private double	DsigmCosAng2Dx2;
	private double	DsigmCosAng2Dy2;
	private double	DsigmCosAng2Dz2;
	private double	DsigmCosAng2Dx3;
	private double	DsigmCosAng2Dy3;
	private double	DsigmCosAng2Dz3;
	private double	DsigmCosAng2Dx4;
	private double	DsigmCosAng2Dy4;
	private double	DsigmCosAng2Dz4;

/**		
 * The constructor parameters:
 * dm - The distance matrix used in this object. All the atoms given as parameters to the 'updateAndEvaluateAtoms1234' 
 * method must be represented in this object.
 * sigmoidBeginsWithH,sigmoidEndsWithH - The transition angles (in degrees) for the HB sigmoid when the hydrogen 
 * atom is known explicitly. Above sigmoidEndsWithH the sigmoid is given a value of 1.0 . Bellow sigmoidBeginsWithH
 * the sigmoid is given a value of 0.0 . In between it raises smoothly by cubic spline.
 * sigmoidBeginsNoH,sigmoidEndsNoH - The same as above, only for HB angle sigmoids when the hydrogen atom is not 
 * given explicitly.
 **/
	public SolvateHBAngle(DistanceMatrix dm , 
							double sigmoidBeginsWithH,
							double sigmoidEndsWithH,
							double sigmoidBeginsNoH,
							double sigmoidEndsNoH) {
		this.dm = dm;
		this.sigmoidBeginsNoH = Math.cos(sigmoidEndsNoH*Math.PI/180.0); // note the switch end-begin. cos is a monotonicly DECREASING function
		this.sigmoidEndsNoH = Math.cos(sigmoidBeginsNoH*Math.PI/180.0);
		this.sigmoidBeginsWithH = Math.cos(sigmoidEndsWithH*Math.PI/180.0);
		this.sigmoidEndsWithH = Math.cos(sigmoidBeginsWithH*Math.PI/180.0);
	}
	
	public final double	hbAngScore() {return hbAngScore;}
	public final boolean zeroDerivative() {return  zeroDerivative;}
	public final double	DhbAngScoreDx1() {return    DhbAngScoreDx1;}
	public final double	DhbAngScoreDy1() {return    DhbAngScoreDy1;}
	public final double	DhbAngScoreDz1() {return    DhbAngScoreDz1;}
	public final double	DhbAngScoreDx2() {return    DhbAngScoreDx2;}
	public final double	DhbAngScoreDy2() {return    DhbAngScoreDy2;}
	public final double	DhbAngScoreDz2() {return    DhbAngScoreDz2;}
	public final double	DhbAngScoreDx3() {return    DhbAngScoreDx3;}
	public final double	DhbAngScoreDy3() {return    DhbAngScoreDy3;}
	public final double	DhbAngScoreDz3() {return    DhbAngScoreDz3;}
	public final double	DhbAngScoreDx4() {return    DhbAngScoreDx4;}
	public final double	DhbAngScoreDy4() {return    DhbAngScoreDy4;}
	public final double	DhbAngScoreDz4() {return    DhbAngScoreDz4;}

	
	public final void updateAndEvaluateAtoms1234(Atom a1,Atom a2,Atom a3,Atom a4) {	
		if (a2.isHydrogen || a3.isHydrogen) {
			sigmoidBegins = sigmoidBeginsWithH;
			sigmoidEnds = sigmoidEndsWithH;
		}
		else {
			sigmoidBegins = sigmoidBeginsNoH;
			sigmoidEnds = sigmoidEndsNoH;			
		}
		
		updateAndEvaluateAtoms123(a1,a2,a3);
		updateAndEvaluateAtoms234(a2,a3,a4);
		
		hbAngScore = sigmCosAng1 * sigmCosAng2;
		if (derivativeAng1Zero && derivativeAng2Zero) {
			zeroDerivative = true;
			DhbAngScoreDx1 = DhbAngScoreDy1 = DhbAngScoreDz1 =
			DhbAngScoreDx2 = DhbAngScoreDy2 = DhbAngScoreDz2 =
			DhbAngScoreDx3 = DhbAngScoreDy3 = DhbAngScoreDz3 =
			DhbAngScoreDx4 = DhbAngScoreDy4 = DhbAngScoreDz4 = 0.0;
		}
		else {
			zeroDerivative = false;
			DhbAngScoreDx1 = DsigmCosAng1Dx1*sigmCosAng2; 
			DhbAngScoreDy1 = DsigmCosAng1Dy1*sigmCosAng2; 
			DhbAngScoreDz1 = DsigmCosAng1Dz1*sigmCosAng2; 
			DhbAngScoreDx2 = DsigmCosAng1Dx2*sigmCosAng2 + DsigmCosAng2Dx2*sigmCosAng1; 
			DhbAngScoreDy2 = DsigmCosAng1Dy2*sigmCosAng2 + DsigmCosAng2Dy2*sigmCosAng1; 
			DhbAngScoreDz2 = DsigmCosAng1Dz2*sigmCosAng2 + DsigmCosAng2Dz2*sigmCosAng1; 
			DhbAngScoreDx3 = DsigmCosAng1Dx3*sigmCosAng2 + DsigmCosAng2Dx3*sigmCosAng1; 
			DhbAngScoreDy3 = DsigmCosAng1Dy3*sigmCosAng2 + DsigmCosAng2Dy3*sigmCosAng1; 
			DhbAngScoreDz3 = DsigmCosAng1Dz3*sigmCosAng2 + DsigmCosAng2Dz3*sigmCosAng1; 
			DhbAngScoreDx4 = DsigmCosAng2Dx4*sigmCosAng1; 
			DhbAngScoreDy4 = DsigmCosAng2Dy4*sigmCosAng1; 
			DhbAngScoreDz4 = DsigmCosAng2Dz4*sigmCosAng1;
		}
	}
	
	private final void updateAndEvaluateAtoms123(Atom a1,Atom a2,Atom a3) {
		Distance dis1 = dm.distance(a1,a2); // distance on pair (a1-a2)
		Distance dis2 = dm.distance(a3,a2); // distance on pair (a3-a2)	
		
		double cosAng = dis1.dDistanceDx()*dis2.dDistanceDx() + dis1.dDistanceDy()*dis2.dDistanceDy() + dis1.dDistanceDz()*dis2.dDistanceDz();
		
		Sigma.sigma(1.0+cosAng , 2 , 1.0+sigmoidBegins ,1.0+sigmoidEnds, 1.0 , 0.0);
		sigmCosAng1 = Sigma.s;
		if (Sigma.s_tag != 0.0) {
			derivativeAng1Zero = false;
			DsigmCosAng1Dx1 = Sigma.s_tag * dis1.invDistance() * (dis2.dDistanceDx() - cosAng*dis1.dDistanceDx());
			DsigmCosAng1Dy1 = Sigma.s_tag * dis1.invDistance() * (dis2.dDistanceDy() - cosAng*dis1.dDistanceDy());
			DsigmCosAng1Dz1 = Sigma.s_tag * dis1.invDistance() * (dis2.dDistanceDz() - cosAng*dis1.dDistanceDz());
			DsigmCosAng1Dx3 = Sigma.s_tag * dis2.invDistance() * (dis1.dDistanceDx() - cosAng*dis2.dDistanceDx());
			DsigmCosAng1Dy3 = Sigma.s_tag * dis2.invDistance() * (dis1.dDistanceDy() - cosAng*dis2.dDistanceDy());
			DsigmCosAng1Dz3 = Sigma.s_tag * dis2.invDistance() * (dis1.dDistanceDz() - cosAng*dis2.dDistanceDz());
			DsigmCosAng1Dx2 = -DsigmCosAng1Dx1 - DsigmCosAng1Dx3;
			DsigmCosAng1Dy2 = -DsigmCosAng1Dy1 - DsigmCosAng1Dy3;
			DsigmCosAng1Dz2 = -DsigmCosAng1Dz1 - DsigmCosAng1Dz3;
		}
		else {
			derivativeAng1Zero = true;
			DsigmCosAng1Dx1 = DsigmCosAng1Dy1 = DsigmCosAng1Dz1 = DsigmCosAng1Dx3 = DsigmCosAng1Dy3 =
			DsigmCosAng1Dz3 = DsigmCosAng1Dx2 = DsigmCosAng1Dy2 = DsigmCosAng1Dz2 = 0.0;
		}			
	}
	
	private final void updateAndEvaluateAtoms234(Atom a2,Atom a3,Atom a4) {
		Distance dis1 = dm.distance(a2,a3); // distance on pair (a2-a3)
		Distance dis2 = dm.distance(a4,a3); // distance on pair (a4-a3)	
		
		double cosAng = dis1.dDistanceDx()*dis2.dDistanceDx() + dis1.dDistanceDy()*dis2.dDistanceDy() + dis1.dDistanceDz()*dis2.dDistanceDz();
		
		Sigma.sigma(1.0+cosAng , 2 , 1.0+sigmoidBegins ,1.0+sigmoidEnds, 1.0 , 0.0);
		sigmCosAng2 = Sigma.s;
		if (Sigma.s_tag != 0.0) {
			derivativeAng2Zero = false;
			DsigmCosAng2Dx2 = Sigma.s_tag * dis1.invDistance() * (dis2.dDistanceDx() - cosAng*dis1.dDistanceDx());
			DsigmCosAng2Dy2 = Sigma.s_tag * dis1.invDistance() * (dis2.dDistanceDy() - cosAng*dis1.dDistanceDy());
			DsigmCosAng2Dz2 = Sigma.s_tag * dis1.invDistance() * (dis2.dDistanceDz() - cosAng*dis1.dDistanceDz());
			DsigmCosAng2Dx4 = Sigma.s_tag * dis2.invDistance() * (dis1.dDistanceDx() - cosAng*dis2.dDistanceDx());
			DsigmCosAng2Dy4 = Sigma.s_tag * dis2.invDistance() * (dis1.dDistanceDy() - cosAng*dis2.dDistanceDy());
			DsigmCosAng2Dz4 = Sigma.s_tag * dis2.invDistance() * (dis1.dDistanceDz() - cosAng*dis2.dDistanceDz());
			DsigmCosAng2Dx3 = -DsigmCosAng2Dx2 - DsigmCosAng2Dx4;
			DsigmCosAng2Dy3 = -DsigmCosAng2Dy2 - DsigmCosAng2Dy4;
			DsigmCosAng2Dz3 = -DsigmCosAng2Dz2 - DsigmCosAng2Dz4;
		}
		else {
			derivativeAng2Zero = true;
			DsigmCosAng2Dx4 = DsigmCosAng2Dy4 = DsigmCosAng2Dz4 = DsigmCosAng2Dx3 = DsigmCosAng2Dy3 =
			DsigmCosAng2Dz3 = DsigmCosAng2Dx2 = DsigmCosAng2Dy2 = DsigmCosAng2Dz2 = 0.0;
		}			
	}
}	
