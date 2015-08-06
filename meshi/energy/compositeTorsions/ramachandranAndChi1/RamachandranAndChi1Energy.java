package meshi.energy.compositeTorsions.ramachandranAndChi1;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.energy.compositeTorsions.ResidueTorsionsList;
import meshi.geometry.DistanceMatrix;

public class RamachandranAndChi1Energy
	extends SimpleEnergyTerm
	implements CompositeTorsionsDefinitions {

	public RamachandranAndChi1Energy() {}

	public RamachandranAndChi1Energy(
			ResidueTorsionsList residueTorsionsList,
			DistanceMatrix distanceMatrix,
			RamachandranAndChi1ParametersList cppl,
			double weight,
			String comment) {
		super( toArray(distanceMatrix, residueTorsionsList), cppl, weight );
		
		this.comment = comment;
		createElementsList( residueTorsionsList );
	}
	
	public EnergyElement createElement(Object baseElement, Parameters parameters) {
		ResidueTorsions resTorsions =
			(ResidueTorsions) baseElement;
		return new RamachandranAndChi1EnergyElement(
					resTorsions, parameters, weight );
	}
	
}
