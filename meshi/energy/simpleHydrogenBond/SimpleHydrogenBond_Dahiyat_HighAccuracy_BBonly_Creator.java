package meshi.energy.simpleHydrogenBond;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.hydrogenBond.Dahiyat.DahiyatHighAccuracyParamaters;
import meshi.geometry.hydrogenBond.Dahiyat.DahiyatHydrogenBondListBBonly;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;

public class SimpleHydrogenBond_Dahiyat_HighAccuracy_BBonly_Creator extends EnergyCreator {
	
	private boolean toCalcDerivatives = true;
	private DahiyatHighAccuracyParamaters parameters = null;

	public SimpleHydrogenBond_Dahiyat_HighAccuracy_BBonly_Creator(double weight , boolean toCalcDerivatives) {
		super(weight);
		this.toCalcDerivatives = toCalcDerivatives;
		parameters = new DahiyatHighAccuracyParamaters();
	}

	public SimpleHydrogenBond_Dahiyat_HighAccuracy_BBonly_Creator(double weight) {
		this(weight, true);
	}

	public AbstractEnergy createEnergyTerm(Protein protein,
			DistanceMatrix distanceMatrix, CommandList commands) {
		return new SimpleHydrogenBondEnergy(new DahiyatHydrogenBondListBBonly(distanceMatrix, protein.atoms(), parameters),
				distanceMatrix,
				weight(),
				toCalcDerivatives);
	}
	
}
