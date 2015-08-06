package meshi.energy.compositeTorsions.ramachandranAndChi1;

import meshi.energy.Parameters;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.energy.compositeTorsions.SplinedPolynomial;
import meshi.energy.compositeTorsions.SplinedPolynomialsLoader;
import meshi.parameters.Residues;

public class RamachandranAndChi13DParameters extends Parameters
	implements CompositeTorsionsDefinitions, Residues {

	/** residue type */
	private final int residueType;

	/** propensity polynomial for this residue */
	private SplinedPolynomial propPolynomial;
	
	public RamachandranAndChi13DParameters( int residueType,
			SplinedPolynomialsLoader spl ) {
		this.residueType = residueType;

		propPolynomial = spl.findPolynomial( residueType, POLYNOMIAL_PHI_PSI_CHI_1_TORSIONS, ALL );
	}

	/** returns rotamer's residue type */
	public int getResidueType() {
		return residueType;
	}
	
	
	public double evaluate( int derivVar, ResidueTorsions resTorsions ) {
		double phi   = resTorsions.getTorsion(PHI).torsion();
		double psi   = resTorsions.getTorsion(PSI).torsion();
		double chi_1 = resTorsions.getTorsion(CHI_1).torsion();
		
		switch (derivVar) {
			case 0: break;
			case 1: derivVar = PHI; break;
			case 2: derivVar = PSI; break;
			case 3: derivVar = CHI_1; break;
		};

		return propPolynomial.value( derivVar, phi, psi, chi_1 );
	}
}
