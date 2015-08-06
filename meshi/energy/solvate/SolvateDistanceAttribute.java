package meshi.energy.solvate;
import meshi.util.MeshiAttribute;

/** 
 * In the SolvateEnergy evaluation function these fields are to be recalculated for EVERY 
 * distance in the non-bonded list. They are used several times in the solvate evaluation
 * process, and we would like to calculate them once. To this end we attach this special 
 * class as an attribute on each Distance instance. 
 * For each distance we calculate two sigmoid values, and an hydrogen-bond strength value: 
 * sigmCa1 - The carbon index of atom 2 on atom 1. If atom 2 is not a carbon then this value
 *           should be zero. If atom 2 is a carbon then this value should be ~1.0 if atom 2 
 *           is spatially near atom 1. This index drops sigmoidally to zero the farther 
 *           atom 2 is.  
 * sigmCa2 - The same as sigmCa1, except detailing the affect of atom 1 on atom 2.
 * sigmHBa1 - The hydrogen bond (HB) strength between atom 1 and 2. If atom 2 can not form HB 
 *           with atom 1, because of its chemical type then this value should be zero.
 *           If atom 2 can form HB with atom 1 then this value should be ~1.0 if atom 1 and 2
 *           are sufficiently close to create a HB, and if their orientation (defined
 *           also by their base atoms - see DahiyatHydrogenBond) permit hydrogen bonding.
 *           This value drops steeply to zero if the conditions to hydrogen bonding 
 *           are violated.
 *              
 * Also provided are the carbon sigmoid values derivative relatives to the atom coordinates. They 
 * have the general form: dsigmCa{1/2}d{x/y/z}{1/2}
 *   
 **/


public class SolvateDistanceAttribute  implements MeshiAttribute {

     public SolvateDistanceAttribute() {}

     public int key() {return SOLVATE_ALL_ATOM_ATTRIBUTE;}

	 // Hydogen involvement
     public boolean isDisInvolvesHydrogens=false;

	 // Fields relevent to polar pairs 
     public boolean isDisBetween2Polars=false;   // Is this a polar pair?	
     public double sigmHB;  // The hydrogen bond value (including distance and angular dependence
     public double saltBridgeFactorA1 = 1.0;     // Salt bridges vs. hydrogen bond strengths
     public double saltBridgeFactorA2 = 1.0;     
     public double saltBridgeFactorForHBenergyA1 = 1.0;     
     public double saltBridgeFactorForHBenergyA2 = 1.0;     


     // Carbon sigmoids and related derivatives
     public double sigmCa1;
     public double dsigmCa1dx1;
     public double dsigmCa1dy1;
     public double dsigmCa1dz1;
     public double dsigmCa1dx2;
     public double dsigmCa1dy2;
     public double dsigmCa1dz2;
     public double sigmCa2;
     public double dsigmCa2dx1;
     public double dsigmCa2dy1;
     public double dsigmCa2dz1;
     public double dsigmCa2dx2;
     public double dsigmCa2dy2;
     public double dsigmCa2dz2;



	public final void resetSigmVals() {
		if (isDisBetween2Polars)
			sigmHB = 0.0;
		else { 
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
     	}
	}

}