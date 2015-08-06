package meshi.energy.compositeTorsions.smoothRotamerLibrary;

import meshi.energy.Parameters;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.energy.compositeTorsions.SplinedPolynomial;
import meshi.parameters.Residues;

/** Parameters used with the SmoothRotamerLibraryEnergy term. Each parameters
 * object contains all the information needed to calculate the term's energy
 * value for a given residue type, between none and six polynomials.
 * Because of the complexity of the energy term, there are several classes
 * that inherit from this parent class, one for each number of possible
 * sidechain torsion angles.
 * 
 * @author El-ad David Amir
 *
 */
public abstract class SmoothRotamerLibraryParameters extends Parameters 
	implements CompositeTorsionsDefinitions, Residues {

	/** residue type */
	private final int residueType;
	
	/** polynomials used */
	protected SplinedPolynomial polynomials[];
	
	/* The following types of polynomials are used:
	 * no. of chi
	 * torsion
	 * angles      polynomials used
	 * 0           none (energy function irrelevant for Glycine and Alanine)
	 * 1           POLYNOMIAL_PHI_PSI, _PHI_PSI_CHI_1
	 * 2           _PHI_PSI, _PHI_PSI_CHI_1, _CHI_1_CHI_2, _CHI_1
	 * 3           _PHI_PSI, _PHI_PSI_CHI_1, _CHI_1_CHI_2, _CHI_1_CHI_3, _CHI_1
	 * 4           _PHI_PSI, _PHI_PSI_CHI_1, _CHI_1_CHI_2, _CHI_1_CHI_3,
	 *             _CHI_1_CHI_4, _CHI_1
	 */
	
	/** Updates residue type of parameters.
	 */
	public SmoothRotamerLibraryParameters( int residueType ) {
		this.residueType = residueType;
		
		if( !legalResidueType() )
			throw new RuntimeException( "parameter residue type mismatch" );
	}
	
	/** returns rotamer's residue type */
	public int getResidueType() {
		return residueType;
	}
	
	/** verifies residue type is a legal residue types for class. */
	protected abstract boolean legalResidueType();
	
	/** retrieves specific polynomial. */
	public SplinedPolynomial getSplinedPolynomial( int spType ) {
		return polynomials[ spType ];
	}
	
	/** evaluates parameters for residue torsions. That is, calculate the
	 * energy value and derivatives. While the parameters class in MESHI
	 * is usually just a static storage, the complex nature of the energy
	 * calculations with the smooth rotamer library allows easier
	 * calculations here.
	 * Note that throughout evaluate() it is assumed that all the torsion
	 * angles required for calculation exist. This is most crucial for
	 * the first and last residue of the protein, which may not have phi
	 * and psi, respectively.
	 * 
	 * @param derivVar torsion angle to be derived in calculation, zero for none.
	 * @param resTorsions all of the calculated residue's torsion angles
	 */
	public abstract double evaluate( int derivVar, ResidueTorsions resTorsions );
	
	/** Reports contents of all polynomials. */
	public void reportPolynomialsList() {
		if( polynomials == null )
			return;
		
		for( SplinedPolynomial splinedPolynomial : polynomials )
			System.out.println( splinedPolynomial );
	}
	
	/** Converts attributes of SmoothRotamerLibraryParameters to string. */
	public String toString() {
		return "residueType: " + residueType + " num. of polynomials: " +
				(polynomials==null ? 0 : polynomials.length);
	}
}
