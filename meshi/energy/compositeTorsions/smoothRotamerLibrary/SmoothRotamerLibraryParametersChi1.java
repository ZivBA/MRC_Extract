package meshi.energy.compositeTorsions.smoothRotamerLibrary;

import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.energy.compositeTorsions.SplinedPolynomial;
import meshi.energy.compositeTorsions.SplinedPolynomialsLoader;

/** Smooth rotamer library parameters for amino acids with one
 * sidechain torsion angle: CYS, SER, THR, VAL.
 * 
 * @author El-ad David Amir
 *
 */
public class SmoothRotamerLibraryParametersChi1 extends
		SmoothRotamerLibraryParameters {
	
	public SmoothRotamerLibraryParametersChi1(
			int residueType, SplinedPolynomialsLoader spl ) {
		super( residueType );

		/* assemble polynomials for amino acid with one sidechain
		 * torsion angle
		 */
		polynomials = new SplinedPolynomial[2];
		polynomials[POLYNOMIAL_PHI_PSI] =
			spl.findPolynomial( residueType, POLYNOMIAL_PHI_PSI_TORSIONS, ALL );
		polynomials[POLYNOMIAL_PHI_PSI_CHI_1] =
			spl.findPolynomial( residueType, POLYNOMIAL_PHI_PSI_CHI_1_TORSIONS, ALL );
	}

	protected boolean legalResidueType() {
		return (getResidueType() == CYS || getResidueType() == SER ||
				getResidueType() == THR || getResidueType() == VAL);
	}
	
	public double evaluate(int derivVar, ResidueTorsions resTorsions) {
		/* evaluate energy of parameter. generally, can be phrased as:
		 * Pr(rot) = Pr(chi_1|phi,psi)
		 *         = Pr(phi,psi,chi_1) / Pr(phi,psi)
		 *         --------->>>
		 * En(rot) = En(phi,psi,chi_1) - En(phi,psi)
		 */
		double phi   = resTorsions.getTorsion(PHI).torsion();
		double psi   = resTorsions.getTorsion(PSI).torsion();
		double chi_1 = resTorsions.getTorsion(CHI_1).torsion();
		SplinedPolynomial enpp1 =
			polynomials[POLYNOMIAL_PHI_PSI_CHI_1];
		SplinedPolynomial enpp =
			polynomials[POLYNOMIAL_PHI_PSI];
		
		if( derivVar == 0 )
			return
				enpp1.value( 0, phi, psi, chi_1 ) -
				enpp.value( 0, phi, psi );
		else if( derivVar == PHI )
			return
				enpp1.value( 1, phi, psi, chi_1 ) -
				enpp.value( 1, phi, psi );
		else if( derivVar == PSI )
			return
				enpp1.value( 2, phi, psi, chi_1 ) -
				enpp.value( 2, phi, psi );
		else if( derivVar == CHI_1 )
			return
				enpp1.value( 3, phi, psi, chi_1 );
		else
			throw new RuntimeException( "derived torsion angle identifier doesn't exist in these parameters");
	}

}
