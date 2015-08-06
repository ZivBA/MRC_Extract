package meshi.energy.simpleHPterm;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;

public class SimpleHPCreator extends EnergyCreator {

	private double weightHydrophobic;
	private double weightHydrophilic;
	private double radiusHydrophobic;
	private double radiusHydrophilic;
	private boolean toCalcDerivatives;
	
	public SimpleHPCreator(double weightHydrophobic, double weightHydrophilic, double radiusHydrophobic, double radiusHydrophilic, boolean toCalcDerivatives) {
		super(weightHydrophobic);
		this.weightHydrophobic = weightHydrophobic;
		this.weightHydrophilic = weightHydrophilic;
		this.radiusHydrophobic = radiusHydrophobic;
		this.radiusHydrophilic = radiusHydrophilic;
		this.toCalcDerivatives = toCalcDerivatives;
	}
	
	public void setHydrophobicWeight(double newWeight) {
		weightHydrophobic = newWeight;
	}
	public void setHydrophilicWeight(double newWeight) {
		weightHydrophilic = newWeight;
	}

	public AbstractEnergy createEnergyTerm(Protein protein,
			DistanceMatrix distanceMatrix, CommandList commands) {
		return new SimpleHP(distanceMatrix,
				weightHydrophobic,
				weightHydrophilic,
				radiusHydrophobic,
				radiusHydrophilic,
				toCalcDerivatives);
	}
	

}
