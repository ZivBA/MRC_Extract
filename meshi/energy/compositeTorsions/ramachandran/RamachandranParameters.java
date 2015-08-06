package meshi.energy.compositeTorsions.ramachandran;

import meshi.energy.Parameters;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.energy.compositeTorsions.SplinedPolynomial;
import meshi.energy.compositeTorsions.SplinedPolynomialsLoader;
import meshi.parameters.Residues;

public class RamachandranParameters extends Parameters
	implements CompositeTorsionsDefinitions, Residues {

	/** residue type */
	private final int residueType;

	/** polynomial for this residue */
	private SplinedPolynomial propPolynomial;
	
	public RamachandranParameters( int residueType,
			SplinedPolynomialsLoader spl ) {
		this.residueType = residueType;

		propPolynomial = spl.findPolynomial( residueType, POLYNOMIAL_PHI_PSI_TORSIONS, ALL );
	}

	/** returns rotamer's residue type */
	public int getResidueType() {
		return residueType;
	}
	
	
	public double evaluate( int derivVar, ResidueTorsions resTorsions ) {
		double phi   = resTorsions.getTorsion(PHI).torsion();
		double psi   = resTorsions.getTorsion(PSI).torsion();
		
		switch (derivVar) {
			case 0: break;
			case 1: derivVar = PHI; break;
			case 2: derivVar = PSI; break;
		};
		
		return propPolynomial.value( derivVar, phi, psi);
	}
}
