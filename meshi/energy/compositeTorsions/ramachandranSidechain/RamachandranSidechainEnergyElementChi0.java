package meshi.energy.compositeTorsions.ramachandranSidechain;

import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;

/** Ramachandran/Sidechain energy element for amino acids with no
 * sidechain torsion angles: ALA, GLY.
 * 
 * @author El-ad David Amir
 *
 */
public class RamachandranSidechainEnergyElementChi0
	extends RamachandranSidechainEnergyElement
	implements CompositeTorsionsDefinitions {

	public RamachandranSidechainEnergyElementChi0(
			ResidueTorsions residueTorsions,
			RamachandranSidechainParameters rsp,
			double weight ) {
		super( residueTorsions, rsp, weight );
	}
	
	protected boolean legalResidueType() {
		return (residueTorsions.getResidueType() == ALA ||
				residueTorsions.getResidueType() == GLY);
	}

	public double evaluate() {
		/* verify energy element is not frozen */
		if( frozen() ) return 0.0;
		
		/* calcualte energy and derivative */
		double energy       = rsp.evaluate( 0, residueTorsions );
		double phi_deriv	= rsp.evaluate( PHI, residueTorsions );
		double psi_deriv	= rsp.evaluate( PSI, residueTorsions );
		
		/* apply weight */
		energy *= weight;
		phi_deriv *= weight;
		psi_deriv *= weight;

		/* apply force to torsions */
		residueTorsions.applyForce( PHI, -phi_deriv );
		residueTorsions.applyForce( PSI, -psi_deriv );

		monitor( energy );
		return energy;		
	}
}
