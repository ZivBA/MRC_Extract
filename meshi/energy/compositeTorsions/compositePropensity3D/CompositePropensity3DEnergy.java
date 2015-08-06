package meshi.energy.compositeTorsions.compositePropensity3D;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.energy.compositeTorsions.ResidueTorsionsList;
import meshi.geometry.DistanceMatrix;

public class CompositePropensity3DEnergy
	extends SimpleEnergyTerm
	implements CompositeTorsionsDefinitions {

	public CompositePropensity3DEnergy() {}

	public CompositePropensity3DEnergy(
			ResidueTorsionsList residueTorsionsList,
			DistanceMatrix distanceMatrix,
			CompositePropensity3DParametersList cppl,
			double weight,
			String comment) {
		super( toArray(distanceMatrix, residueTorsionsList), cppl, weight );
		
		this.comment = comment;
		createElementsList( residueTorsionsList );
	}
	
	public EnergyElement createElement(Object baseElement, Parameters parameters) {
		ResidueTorsions resTorsions =
			(ResidueTorsions) baseElement;
		return new CompositePropensity3DEnergyElement(
					resTorsions, parameters, weight );
	}
	
}
