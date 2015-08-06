package meshi.energy.compositeTorsions.ramachandran;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.energy.compositeTorsions.ResidueTorsionsList;
import meshi.geometry.DistanceMatrix;

public class RamachandranEnergy
	extends SimpleEnergyTerm
	implements CompositeTorsionsDefinitions {

	public RamachandranEnergy() {}

	public RamachandranEnergy(
			ResidueTorsionsList residueTorsionsList,
			DistanceMatrix distanceMatrix,
			RamachandranParametersList cppl,
			double weight,
			String comment) {
		super( toArray(distanceMatrix, residueTorsionsList), cppl, weight );
		
		this.comment = comment;
		createElementsList( residueTorsionsList );
	}
	
	public EnergyElement createElement(Object baseElement, Parameters parameters) {
		ResidueTorsions resTorsions =
			(ResidueTorsions) baseElement;
		RamachandranParameters cpp =
			(RamachandranParameters) parameters;
		
		return new RamachandranEnergyElement(
					resTorsions, cpp, weight );
	}
	
}
