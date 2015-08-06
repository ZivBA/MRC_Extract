package meshi.energy.compositeTorsions.smoothRotamerLibrary;

import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;

/** Smooth rotamer library energy element for amino acids with three
 * sidechain torsion angles: GLU, MET, GLN.
 * 
 * @author El-ad David Amir
 *
 */
public class SmoothRotamerLibraryEnergyElementChi3
	extends SmoothRotamerLibraryEnergyElement
	implements CompositeTorsionsDefinitions {

	public SmoothRotamerLibraryEnergyElementChi3(
			ResidueTorsions residueTorsions,
			SmoothRotamerLibraryParameters srlp,
			double weight ) {
		super( residueTorsions, srlp, weight );
	}

	protected boolean legalResidueType() {
		return (residueTorsions.getResidueType() == GLU ||
				residueTorsions.getResidueType() == MET ||
				residueTorsions.getResidueType() == GLN);
	}

	public double evaluate() {
		/* verify energy element is not frozen */
		if( frozen() ) return 0.0;
		
		/* calcualte energy and derivative */
		double energy       = srlp.evaluate( 0, residueTorsions );
		double phi_deriv	= srlp.evaluate( PHI, residueTorsions );
		double psi_deriv	= srlp.evaluate( PSI, residueTorsions );
		double chi_1_deriv  = srlp.evaluate( CHI_1, residueTorsions );
		double chi_2_deriv  = srlp.evaluate( CHI_2, residueTorsions );
		double chi_3_deriv  = srlp.evaluate( CHI_3, residueTorsions );
		
		/* apply weight */
		energy *= weight;
		phi_deriv *= weight;
		psi_deriv *= weight;
		chi_1_deriv *= weight;
		chi_2_deriv *= weight;
		chi_3_deriv *= weight;

		/* apply force to torsions */
		residueTorsions.applyForce( PHI, -phi_deriv );
		residueTorsions.applyForce( PSI, -psi_deriv );
		residueTorsions.applyForce( CHI_1, -chi_1_deriv );
		residueTorsions.applyForce( CHI_2, -chi_2_deriv );
		residueTorsions.applyForce( CHI_3, -chi_3_deriv );

		monitor( energy, chi_1_deriv, chi_2_deriv, chi_3_deriv );
		return energy;
	}

}
