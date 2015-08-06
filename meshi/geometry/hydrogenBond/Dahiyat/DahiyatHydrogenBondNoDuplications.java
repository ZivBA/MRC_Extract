package meshi.geometry.hydrogenBond.Dahiyat;

import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.AtomList;
import meshi.parameters.AtomTypes;
import meshi.util.mathTools.HB_Sigma;

/**
 * This hydrogen bond is to be used with the 'DahiyatHydrogenBondListNoDuplications' class. See the comments there for more information.
 * 
 * @author Nir
 *
 */

public class DahiyatHydrogenBondNoDuplications extends DahiyatHydrogenBond implements AtomTypes {

	// Is this hydrogen bond involving ASP, GLU or ARG functional groups?
	private boolean isPolar1Relevant = false;
	private boolean isPolar2Relevant = false;

	// This weight can be zeroed in case of a duplication, for the weaker bond.
	private double duplicationWeight = 1.0;
	
	public DahiyatHydrogenBondNoDuplications() {
    	throw new RuntimeException("\nERROR: without parameters the hydrogen bonds cannot be formed.\n");
	}

	public DahiyatHydrogenBondNoDuplications(DistanceMatrix dm,
			AtomList atomList, HB_Sigma sigmaDis,
			DahiyatParametersInterface angleParameters) {
		super(dm, atomList, sigmaDis, angleParameters);
		if ((getFirstPolar().type==DOD) || (getFirstPolar().type==EOE) || (getFirstPolar().type==RNH))
			isPolar1Relevant = true;
		if ((getSecondPolar().type==DOD) || (getSecondPolar().type==EOE) || (getSecondPolar().type==RNH))
			isPolar2Relevant = true;
	}
	
	public double hbVal() {
		return super.hbVal()* duplicationWeight;
	}

	public double hbValWithoutDuplicationWeight() {
		return super.hbVal();
	}
	
	protected void setDuplicationWeight(double setTo) {
		duplicationWeight = setTo;
	}

	protected double getDuplicationWeight() {
		return duplicationWeight;
	}

	public void applyForcesToAtoms(double factor) {
		super.applyForcesToAtoms(factor*duplicationWeight);
	}
	
	public boolean isPolar1Relevant() { 
		return isPolar1Relevant;
	}

	public boolean isPolar2Relevant() { 
		return isPolar2Relevant;
	}
}
