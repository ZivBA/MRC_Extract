package meshi.energy.compositeTorsions.compositePropensity3D;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.parameters.Residues;

public class CompositePropensity3DEnergyElement
	extends EnergyElement
	implements Residues, CompositeTorsionsDefinitions {

	public final ResidueTorsions residueTorsions;
	public final CompositePropensity2DParameters cpp2D;
	public final CompositePropensity3DParameters cpp3D;
	protected double weight;
	public final boolean is3D;

	public CompositePropensity3DEnergyElement(
			ResidueTorsions residueTorsions,
			Parameters cpp,
			double weight ) {
		this.residueTorsions = residueTorsions;
		if (residueTorsions.getTorsion( CHI_1 ) != null) {
			cpp2D=null;
			cpp3D = (CompositePropensity3DParameters) cpp;
			is3D = true;
		}
		else {
			cpp2D = (CompositePropensity2DParameters) cpp;
			cpp3D = null;
			is3D = false;			
		}
		this.weight = weight;
		
		setAtoms();
		updateFrozen();
	}

	protected void setAtoms() {
		int[] interestingTorsions = {PHI,PSI,CHI_1};
		atoms = residueTorsions.getAtoms(interestingTorsions);		
	}

	public double evaluate() {
		double energy = 0.0;

 		/* verify energy element is not frozen */
		if( frozen() ) return energy;
		
		if (is3D) {
		/* calcualte energy and derivative */
		energy       = cpp3D.evaluate( 0, residueTorsions );
		double phi_deriv	= cpp3D.evaluate( PHI, residueTorsions );
		double psi_deriv	= cpp3D.evaluate( PSI, residueTorsions );
		double chi_1_deriv  = cpp3D.evaluate( CHI_1, residueTorsions );
		
		/* apply weight */
		energy *= weight;
		phi_deriv *= weight;
		psi_deriv *= weight;
		chi_1_deriv *= weight;

		/* apply force to torsions */
		residueTorsions.applyForce( PHI, -phi_deriv );
		residueTorsions.applyForce( PSI, -psi_deriv );
		residueTorsions.applyForce( CHI_1, -psi_deriv );
		}
		else {
		/* calcualte energy and derivative */
		energy       = cpp2D.evaluate( 0, residueTorsions );
		double phi_deriv	= cpp2D.evaluate( PHI, residueTorsions );
		double psi_deriv	= cpp2D.evaluate( PSI, residueTorsions );
		
		/* apply weight */
		energy *= weight;
		phi_deriv *= weight;
		psi_deriv *= weight;

		/* apply force to torsions */
		residueTorsions.applyForce( PHI, -phi_deriv );
		residueTorsions.applyForce( PSI, -psi_deriv );
		}

		return energy;		
	}	
}
