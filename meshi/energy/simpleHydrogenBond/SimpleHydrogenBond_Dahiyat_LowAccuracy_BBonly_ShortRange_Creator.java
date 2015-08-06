package meshi.energy.simpleHydrogenBond;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.hydrogenBond.Dahiyat.DahiyatHydrogenBondListBBonlyShortRange;
import meshi.geometry.hydrogenBond.Dahiyat.DahiyatLowAccuracyParameters;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;

public class SimpleHydrogenBond_Dahiyat_LowAccuracy_BBonly_ShortRange_Creator extends EnergyCreator {

	private boolean toCalcDerivatives = true;
	private DahiyatLowAccuracyParameters parameters = null;

	public SimpleHydrogenBond_Dahiyat_LowAccuracy_BBonly_ShortRange_Creator(double weight , boolean toCalcDerivatives) {
		super(weight);
		this.toCalcDerivatives = toCalcDerivatives;
		parameters = new DahiyatLowAccuracyParameters();
	}

	public SimpleHydrogenBond_Dahiyat_LowAccuracy_BBonly_ShortRange_Creator(double weight) {
		this(weight, true);
	}

	public AbstractEnergy createEnergyTerm(Protein protein,
			DistanceMatrix distanceMatrix, CommandList commands) {
		return new SimpleHydrogenBondEnergy(new DahiyatHydrogenBondListBBonlyShortRange(distanceMatrix, protein.atoms(), parameters),
				distanceMatrix,
				weight(),
				toCalcDerivatives);
	}

}
