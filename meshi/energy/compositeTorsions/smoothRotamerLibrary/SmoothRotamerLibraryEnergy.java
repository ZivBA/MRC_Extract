package meshi.energy.compositeTorsions.smoothRotamerLibrary;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.energy.compositeTorsions.ResidueTorsionsList;
import meshi.geometry.DistanceMatrix;

/** Smooth rotamer library energy term. Statistical analysis of
 * residue backbone and sidechain torsions in a large database of
 * residue observations has been smoothed using polynomial spline
 * interpolation.
 * For a given residue the energy value approximates the percentage
 * of finding its current sidechain torsion angles.
 * 
 * @author El-ad David Amir
 *
 */
public class SmoothRotamerLibraryEnergy
	extends SimpleEnergyTerm
	implements CompositeTorsionsDefinitions {

	public SmoothRotamerLibraryEnergy(){}

	public SmoothRotamerLibraryEnergy(
			ResidueTorsionsList residueTorsionsList,
			DistanceMatrix distanceMatrix,
			SmoothRotamerLibraryParametersList srlpl,
			double weight,
			String comment) {
		super( toArray(distanceMatrix, residueTorsionsList), srlpl, weight );
		
		this.comment = comment;
		createElementsList( residueTorsionsList );
	}
	
	public EnergyElement createElement(Object baseElement, Parameters parameters) {
		ResidueTorsions resTorsions =
			(ResidueTorsions) baseElement;
		SmoothRotamerLibraryParameters srlp =
			(SmoothRotamerLibraryParameters) parameters;
		
		switch( NUM_SIDECHAIN_TORSIONS[resTorsions.getResidueType()] ) {
		case 0:
			return new SmoothRotamerLibraryEnergyElementChi0(
					resTorsions, srlp, weight );
		case 1:
			return new SmoothRotamerLibraryEnergyElementChi1(
					resTorsions, srlp, weight );
		case 2:
			return new SmoothRotamerLibraryEnergyElementChi2(
					resTorsions, srlp, weight );
		case 3:
			return new SmoothRotamerLibraryEnergyElementChi3(
					resTorsions, srlp, weight );
		case 4:
			return new SmoothRotamerLibraryEnergyElementChi4(
					resTorsions, srlp, weight );
		default:
			throw new RuntimeException( "unidentified residue type" );
		}
	}

}
