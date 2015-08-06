package meshi.energy.compositeTorsions.smoothRotamerLibrary;

import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;

/** Smooth rotamer library energy element for amino acids with one
 * sidechain torsion angle: CYS, SER, THR, VAL.
 * 
 * @author El-ad David Amir
 *
 */
public class SmoothRotamerLibraryEnergyElementChi1
	extends SmoothRotamerLibraryEnergyElement
	implements CompositeTorsionsDefinitions {

	public SmoothRotamerLibraryEnergyElementChi1(
			ResidueTorsions residueTorsions,
			SmoothRotamerLibraryParameters srlp,
			double weight ) {
		super( residueTorsions, srlp, weight );
	}

	protected boolean legalResidueType() {
		return (residueTorsions.getResidueType() == CYS ||
				residueTorsions.getResidueType() == SER ||
				residueTorsions.getResidueType() == THR ||
				residueTorsions.getResidueType() == VAL);
	}
	
	public double evaluate() {
		/* verify energy element is not frozen */
		if( frozen() ) return 0.0;
		
		/* calcualte energy and derivative */
		double energy       = srlp.evaluate( 0, residueTorsions );
		double phi_deriv	= srlp.evaluate( PHI, residueTorsions );
		double psi_deriv	= srlp.evaluate( PSI, residueTorsions );
		double chi_1_deriv  = srlp.evaluate( CHI_1, residueTorsions );
		
		/* apply weight */
		energy *= weight;
		phi_deriv *= weight;
		psi_deriv *= weight;
		chi_1_deriv *= weight;
		
		/* apply force to torsions */
		residueTorsions.applyForce( PHI, -phi_deriv );
		residueTorsions.applyForce( PSI, -psi_deriv );
		residueTorsions.applyForce( CHI_1, -chi_1_deriv );

		monitor( energy, chi_1_deriv );
		return energy;
	}

}
