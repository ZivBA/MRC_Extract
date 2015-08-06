package meshi.energy.compositeTorsions.smoothRotamerLibrary;

import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.energy.compositeTorsions.SplinedPolynomial;
import meshi.energy.compositeTorsions.SplinedPolynomialsLoader;

/** Smooth rotamer library parameters for amino acids with two
 * sidechain torsion angles: ASP, PHE, HIS, ILE, LEU, ASN, PRO,
 * TRP, TYR.
 * 
 * @author El-ad David Amir
 *
 */
public class SmoothRotamerLibraryParametersChi2 extends
		SmoothRotamerLibraryParameters {

	public SmoothRotamerLibraryParametersChi2(
			int residueType, SplinedPolynomialsLoader spl ) {
		super( residueType );

		/* assemble polynomials for amino acid with two sidechain
		 * torsion angles
		 */
		polynomials = new SplinedPolynomial[4];
		polynomials[POLYNOMIAL_PHI_PSI] =
			spl.findPolynomial( residueType, POLYNOMIAL_PHI_PSI_TORSIONS, ALL );
		polynomials[POLYNOMIAL_PHI_PSI_CHI_1] =
			spl.findPolynomial( residueType, POLYNOMIAL_PHI_PSI_CHI_1_TORSIONS, ALL );
		polynomials[POLYNOMIAL_CHI_1_CHI_2] =
			spl.findPolynomial( residueType, POLYNOMIAL_CHI_1_CHI_2_TORSIONS, ALL );
		polynomials[POLYNOMIAL_CHI_1] =
			spl.findPolynomial( residueType, POLYNOMIAL_CHI_1_TORSIONS, ALL );
	}

	protected boolean legalResidueType() {
		return (getResidueType() == ASP || getResidueType() == PHE ||
				getResidueType() == HIS || getResidueType() == ILE ||
				getResidueType() == LEU || getResidueType() == ASN ||
				getResidueType() == PRO || getResidueType() == TRP ||
				getResidueType() == TYR);
	}

	public double evaluate(int derivVar, ResidueTorsions resTorsions) {
		/* evaluate energy of parameter. generally, can be phrased as:
		 * Pr(rot) = Pr(chi_1|phi,psi) * Pr(chi_2|chi_1)
		 *         = (Pr(phi,psi,chi_1) / Pr(phi,psi)) *
		 *           (Pr(chi_1,chi_2) / Pr(chi_1))
		 *         --------->>>
		 * En(rot) = En(phi,psi,chi_1) - En(phi,psi) +
		 *           En(chi_1,chi_2) - En(chi_1)
		 */
		double phi   = resTorsions.getTorsion(PHI).torsion();
		double psi   = resTorsions.getTorsion(PSI).torsion();
		double chi_1 = resTorsions.getTorsion(CHI_1).torsion();
		double chi_2 = resTorsions.getTorsion(CHI_2).torsion();
		SplinedPolynomial enpp1 =
			polynomials[POLYNOMIAL_PHI_PSI_CHI_1];
		SplinedPolynomial enpp =
			polynomials[POLYNOMIAL_PHI_PSI];
		SplinedPolynomial en12 =
			polynomials[POLYNOMIAL_CHI_1_CHI_2];
		SplinedPolynomial en1 =
			polynomials[POLYNOMIAL_CHI_1];
		
		if( derivVar == 0 )
			return
				enpp1.value( 0, phi, psi, chi_1 ) -
				enpp.value( 0, phi, psi ) +
				en12.value( 0, chi_1, chi_2 ) -
				en1.value( 0, chi_1 );
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
				enpp1.value( 3, phi, psi, chi_1 ) +
				en12.value( 1, chi_1, chi_2 ) -
				en1.value( 1, chi_1 );
		else if( derivVar == CHI_2 )
			return
				en12.value( 2, chi_1, chi_2 );			
		else
			throw new RuntimeException( "derived torsion angle identifier doesn't exist in these parameters");
	}

}
