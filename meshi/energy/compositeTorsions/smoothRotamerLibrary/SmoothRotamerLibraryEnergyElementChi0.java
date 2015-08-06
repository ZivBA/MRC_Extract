package meshi.energy.compositeTorsions.smoothRotamerLibrary;

import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;

/** Smooth rotamer library energy element for amino acids with no
 * sidechain torsion angles: ALA, GLY.
 * 
 * @author El-ad David Amir
 *
 */
public class SmoothRotamerLibraryEnergyElementChi0
	extends SmoothRotamerLibraryEnergyElement
	implements CompositeTorsionsDefinitions {

	public SmoothRotamerLibraryEnergyElementChi0(
			ResidueTorsions residueTorsions,
			SmoothRotamerLibraryParameters srlp,
			double weight ) {
		super( residueTorsions, srlp, weight );
	}
	
	public double evaluate() {
		/* no sidechain, evaluation is always zero. we allow ourselves
		 * a shortcut here, rather than call the parameter's evaluation
		 * method which we know will return 0.
		 */
		return 0.0;
	}

	protected boolean legalResidueType() {
		return (residueTorsions.getResidueType() == ALA ||
				residueTorsions.getResidueType() == GLY);
	}
}
