package meshi.energy.oldSolvate;
import meshi.util.MeshiAttribute;

/** 
 * In the SolvateEnergy evaluation function these fields are to be recalculated for EVERY 
 * distance in the non-bonded list. They are used several times in the solvate evaluation
 * process, and we would like to calculate them once. To this end we attach this special 
 * class as an attribute on each Distance instance. 
 * For each distance we calculate four sigmoid values: two for atom 1 in the distance (a1), 
 * and two for atom 2 in the distance (a2):
 * sigmCa1 - The carbon index of atom 2 on atom 1. If atom 2 is not a carbon then this value
 *           should be zero. If atom 2 is a carbon then this value should be ~1.0 if atom 2 
 *           is spatially near atom 1. This index drops sigmoidally to zero the farther 
 *           atom 2 is.  
 * sigmHBa1 - The hydrogen bond (HB) index of atom 2 on atom 1. If atom 2 can not form HB 
 *           with atom 1, because of its chemical type then this value should be zero.
 *           If atom 2 can form HB with atom 1 then this value should be ~1.0 if atom 1 and 2
 *           are sufficiently close to create a HB, and if their orientation (defined
 *           also by their base atoms - see SolvateHBAngle) permit hydrogen bonding.
 *           This index drops drop steeply to zero if the conditions to hydrogen bonding 
 *           are violated.
 * sigmCa2 - The same as sigmCa1, except detailing the affect of atom 1 on atom 2.
 * sigmHBa2 - The same as sigmHBa1, except detailing the affect of atom 1 on atom 2.
 *   
 * Also provided are the sigmoid values derivative relatives to the atom coordinates. They 
 * have the general form: dsigm{HB/C}a{1/2}d{x/y/z}{1/2/3/4}
 * Because of their orientation dependency, the HB sigmoids involve 4 atoms (the two atoms of 
 * the distance and their respective base atoms). Therefore, the derivatives of the sigmoid 
 * also include atom 3 (the base atom of atom 1), and atom 4 (the base atom of atom 2).  
 *   
 **/


public class SolvateDistanceAttribute  implements MeshiAttribute {

     public SolvateDistanceAttribute() {}

     public int key() {return SOLVATE_ALL_ATOM_ATTRIBUTE;}

     public double sigmCa1;
     public double sigmHBa1;
     public double sigmPHBa1;
     public double sigmCa2;
     public double sigmHBa2;
     public double sigmPHBa2;
     public double dsigmCa1dx1;
     public double dsigmCa1dy1;
     public double dsigmCa1dz1;
     public double dsigmCa1dx2;
     public double dsigmCa1dy2;
     public double dsigmCa1dz2;
     public double dsigmHBa1dx1;
     public double dsigmHBa1dy1;
     public double dsigmHBa1dz1;
     public double dsigmHBa1dx2;
     public double dsigmHBa1dy2;
     public double dsigmHBa1dz2;
     public double dsigmPHBa1dx1;
     public double dsigmPHBa1dy1;
     public double dsigmPHBa1dz1;
     public double dsigmPHBa1dx2;
     public double dsigmPHBa1dy2;
     public double dsigmPHBa1dz2;
     public double dsigmCa2dx1;
     public double dsigmCa2dy1;
     public double dsigmCa2dz1;
     public double dsigmCa2dx2;
     public double dsigmCa2dy2;
     public double dsigmCa2dz2;
     public double dsigmHBa2dx1;
     public double dsigmHBa2dy1;
     public double dsigmHBa2dz1;
     public double dsigmHBa2dx2;
     public double dsigmHBa2dy2;
     public double dsigmHBa2dz2;
     public double dsigmPHBa2dx1;
     public double dsigmPHBa2dy1;
     public double dsigmPHBa2dz1;
     public double dsigmPHBa2dx2;
     public double dsigmPHBa2dy2;
     public double dsigmPHBa2dz2;
     // values arising from the fact that the angles of a HB involve 4 atoms 
     public boolean hbAngleDerivativeNonZero;
     public double dsigmHBa1dx3;
     public double dsigmHBa1dy3;
     public double dsigmHBa1dz3;
     public double dsigmHBa1dx4;
     public double dsigmHBa1dy4;
     public double dsigmHBa1dz4;
     public double dsigmHBa2dx3;
     public double dsigmHBa2dy3;
     public double dsigmHBa2dz3;
     public double dsigmHBa2dx4;
     public double dsigmHBa2dy4;
     public double dsigmHBa2dz4;
     public double dsigmPHBa1dx3;
     public double dsigmPHBa1dy3;
     public double dsigmPHBa1dz3;
     public double dsigmPHBa1dx4;
     public double dsigmPHBa1dy4;
     public double dsigmPHBa1dz4;
     public double dsigmPHBa2dx3;
     public double dsigmPHBa2dy3;
     public double dsigmPHBa2dz3;
     public double dsigmPHBa2dx4;
     public double dsigmPHBa2dy4;
     public double dsigmPHBa2dz4;

	private final void resetHBSigmVals() {
   		sigmHBa1 = 0.0;
   		dsigmHBa1dx1 = 0.0;
   		dsigmHBa1dy1 = 0.0;
   		dsigmHBa1dz1 = 0.0;
   		dsigmHBa1dx2 = 0.0;
   		dsigmHBa1dy2 = 0.0;
   		dsigmHBa1dz2 = 0.0;
   		sigmHBa2 = 0.0;
   		dsigmHBa2dx1 = 0.0;
   		dsigmHBa2dy1 = 0.0;
   		dsigmHBa2dz1 = 0.0;
   		dsigmHBa2dx2 = 0.0;
   		dsigmHBa2dy2 = 0.0;
   		dsigmHBa2dz2 = 0.0;
   		dsigmHBa1dx3 = 0.0;
    	dsigmHBa1dy3 = 0.0;
    	dsigmHBa1dz3 = 0.0;
    	dsigmHBa1dx4 = 0.0;
    	dsigmHBa1dy4 = 0.0;
     	dsigmHBa1dz4 = 0.0;
     	dsigmHBa2dx3 = 0.0;
     	dsigmHBa2dy3 = 0.0;
     	dsigmHBa2dz3 = 0.0;
     	dsigmHBa2dx4 = 0.0;
     	dsigmHBa2dy4 = 0.0;
     	dsigmHBa2dz4 = 0.0;
   		sigmPHBa1 = 0.0;
   		dsigmPHBa1dx1 = 0.0;
   		dsigmPHBa1dy1 = 0.0;
   		dsigmPHBa1dz1 = 0.0;
   		dsigmPHBa1dx2 = 0.0;
   		dsigmPHBa1dy2 = 0.0;
   		dsigmPHBa1dz2 = 0.0;
   		sigmPHBa2 = 0.0;
   		dsigmPHBa2dx1 = 0.0;
   		dsigmPHBa2dy1 = 0.0;
   		dsigmPHBa2dz1 = 0.0;
   		dsigmPHBa2dx2 = 0.0;
   		dsigmPHBa2dy2 = 0.0;
   		dsigmPHBa2dz2 = 0.0;
   		dsigmPHBa1dx3 = 0.0;
    	dsigmPHBa1dy3 = 0.0;
    	dsigmPHBa1dz3 = 0.0;
    	dsigmPHBa1dx4 = 0.0;
    	dsigmPHBa1dy4 = 0.0;
     	dsigmPHBa1dz4 = 0.0;
     	dsigmPHBa2dx3 = 0.0;
     	dsigmPHBa2dy3 = 0.0;
     	dsigmPHBa2dz3 = 0.0;
     	dsigmPHBa2dx4 = 0.0;
     	dsigmPHBa2dy4 = 0.0;
     	dsigmPHBa2dz4 = 0.0;
     	hbAngleDerivativeNonZero = false;
    }


	public final void resetAllSigmVals() {
		sigmCa1 = 0.0;
    	dsigmCa1dx1 = 0.0;
    	dsigmCa1dy1 = 0.0;
    	dsigmCa1dz1 = 0.0;
	    dsigmCa1dx2 = 0.0;
    	dsigmCa1dy2 = 0.0;
    	dsigmCa1dz2 = 0.0;
	    sigmCa2 = 0.0;
	    dsigmCa2dx1 = 0.0;
	    dsigmCa2dy1 = 0.0;
	    dsigmCa2dz1 = 0.0;
	    dsigmCa2dx2 = 0.0;
    	dsigmCa2dy2 = 0.0;
    	dsigmCa2dz2 = 0.0;
		resetHBSigmVals();
	}

	
}