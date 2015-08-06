package meshi.energy.compositeTorsions.ramachandranSidechain;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.energy.compositeTorsions.ResidueTorsionsList;
import meshi.geometry.DistanceMatrix;

/** A Ramachandran plot and sidechain torsions optimization energy term.
 * Statistical analysis of residue backbone and sidechain torsions in a
 * large database of residue observations has been smoothed using
 * polynomial spline interpolation.
 * For a given residue the energy value approximates the percentage
 * of finding its current backbone and sidechain torsion angles.
 * 
 * @author El-ad David Amir
 *
 */
public class RamachandranSidechainEnergy
	extends SimpleEnergyTerm
	implements CompositeTorsionsDefinitions {

	public RamachandranSidechainEnergy() {}

	public RamachandranSidechainEnergy(
			ResidueTorsionsList residueTorsionsList,
			DistanceMatrix distanceMatrix,
			RamachandranSidechainParametersList rspl,
			double weight,
			String comment) {
		super( toArray(distanceMatrix, residueTorsionsList), rspl, weight );
		
		this.comment = comment;
		createElementsList( residueTorsionsList );
	}
	
	public EnergyElement createElement(Object baseElement, Parameters parameters) {
		ResidueTorsions resTorsions =
			(ResidueTorsions) baseElement;
		RamachandranSidechainParameters rsp =
			(RamachandranSidechainParameters) parameters;
		
		switch( NUM_SIDECHAIN_TORSIONS[resTorsions.getResidueType()] ) {
		case 0:
			return new RamachandranSidechainEnergyElementChi0(
					resTorsions, rsp, weight );
		case 1:
			return new RamachandranSidechainEnergyElementChi1(
					resTorsions, rsp, weight );
		case 2:
			return new RamachandranSidechainEnergyElementChi2(
					resTorsions, rsp, weight );
		case 3:
			return new RamachandranSidechainEnergyElementChi3(
					resTorsions, rsp, weight );
		case 4:
			return new RamachandranSidechainEnergyElementChi4(
					resTorsions, rsp, weight );
		default:
			throw new RuntimeException( "unidentified residue type" );
		}
	}
	
}
