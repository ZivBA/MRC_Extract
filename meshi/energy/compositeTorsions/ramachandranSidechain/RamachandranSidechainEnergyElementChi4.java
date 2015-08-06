package meshi.energy.compositeTorsions.ramachandranSidechain;

import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;

/** Ramachandran/Sidechain parameters for amino acids with four
 * sidechain torsion angles: LYS, ARG.
 * 
 * @author El-ad David Amir
 *
 */
public class RamachandranSidechainEnergyElementChi4
	extends RamachandranSidechainEnergyElement
	implements CompositeTorsionsDefinitions {

	public RamachandranSidechainEnergyElementChi4(
			ResidueTorsions residueTorsions,
			RamachandranSidechainParameters rsp,
			double weight ) {
		super( residueTorsions, rsp, weight );
	}
	
	protected boolean legalResidueType() {
		return (residueTorsions.getResidueType() == LYS ||
				residueTorsions.getResidueType() == ARG);
	}

	public double evaluate() {
		/* verify energy element is not frozen */
		if( frozen() ) return 0.0;
		
		/* calcualte energy and derivative */
		double energy       = rsp.evaluate( 0, residueTorsions );
		double phi_deriv	= rsp.evaluate( PHI, residueTorsions );
		double psi_deriv	= rsp.evaluate( PSI, residueTorsions );
		double chi_1_deriv  = rsp.evaluate( CHI_1, residueTorsions );
		double chi_2_deriv  = rsp.evaluate( CHI_2, residueTorsions );
		double chi_3_deriv  = rsp.evaluate( CHI_3, residueTorsions );
		double chi_4_deriv  = rsp.evaluate( CHI_4, residueTorsions );
		
		/* apply weight */
		energy *= weight;
		phi_deriv *= weight;
		psi_deriv *= weight;
		chi_1_deriv *= weight;
		chi_2_deriv *= weight;
		chi_3_deriv *= weight;
		chi_4_deriv *= weight;

		/* apply force to torsions */
		residueTorsions.applyForce( PHI, -phi_deriv );
		residueTorsions.applyForce( PSI, -psi_deriv );
		residueTorsions.applyForce( CHI_1, -chi_1_deriv );
		residueTorsions.applyForce( CHI_2, -chi_2_deriv );
		residueTorsions.applyForce( CHI_3, -chi_3_deriv );
		residueTorsions.applyForce( CHI_4, -chi_4_deriv );

		monitor( energy, chi_1_deriv, chi_2_deriv, chi_3_deriv, chi_4_deriv );
		return energy;
	}

}
