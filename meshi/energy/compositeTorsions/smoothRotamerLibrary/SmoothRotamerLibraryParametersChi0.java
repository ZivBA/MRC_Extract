package meshi.energy.compositeTorsions.smoothRotamerLibrary;

import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.energy.compositeTorsions.SplinedPolynomialsLoader;

/** Smooth rotamer library parameters for amino acids with no
 * sidechain torsion angles: ALA, GLY.
 * 
 * @author El-ad David Amir
 *
 */
public class SmoothRotamerLibraryParametersChi0 extends
		SmoothRotamerLibraryParameters {

	public SmoothRotamerLibraryParametersChi0(
			int residueType, SplinedPolynomialsLoader spl ) {
		super( residueType );
		
		polynomials = null;
	}
	
	protected boolean legalResidueType() {
		return (getResidueType() == ALA || getResidueType() == GLY);
	}
	
	public double evaluate(int derivVar, ResidueTorsions resTorsions) {
		/* no sidechain, evaluation is always zero */
		return 0;
	}

}
