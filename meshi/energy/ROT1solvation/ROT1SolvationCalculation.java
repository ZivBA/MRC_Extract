package meshi.energy.ROT1solvation;

import meshi.energy.ROT1solvation.parameters.AbstractROT1Parameters;

/**
 * This class encapsulate the calculation of the Rot1 solvation.
 * You need to give a reference to a parameter class.
 * 
 * So far I assume that the 'getROT1SolvValueArray' is giving the solvation values for cnc's of:
 * {0,1,2,3,.....}, i.e. equal spacing of 1 for the cnc values.
 * Solvation for (cnc<0) or (cnc>=cncMax) is 0.0
 * 
 * @author Nir
 *
 */
public class ROT1SolvationCalculation {
	
	private double[][] solvArray = null;
	
	public ROT1SolvationCalculation(AbstractROT1Parameters params) {
		solvArray = params.getROT1SolvValueArray();
	}
		
	public double calcRot1(int atomType, double cnc) {
		int ind = (int) cnc;
		if ((ind>=(solvArray[atomType].length-1)) || (ind<0))
			return 0.0;
		return solvArray[atomType][ind]+(cnc-ind)*(solvArray[atomType][ind+1]-solvArray[atomType][ind]);
	}

	public double calcRot1Deriv(int atomType, double cnc) {
		int ind = (int) cnc;
		if ((ind>=(solvArray[atomType].length-1)) || (ind<0))
			return 0.0;
		return (solvArray[atomType][ind+1]-solvArray[atomType][ind]);
	}
		
}
