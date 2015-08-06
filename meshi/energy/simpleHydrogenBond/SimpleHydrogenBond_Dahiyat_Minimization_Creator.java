package meshi.energy.simpleHydrogenBond;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.hydrogenBond.Dahiyat.DahiyatHydrogenBondList;
import meshi.geometry.hydrogenBond.Dahiyat.DahiyatMinimizationParameters;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;

public class SimpleHydrogenBond_Dahiyat_Minimization_Creator extends EnergyCreator {

	private boolean toCalcDerivatives = true;
	private DahiyatMinimizationParameters parameters = null;

	public SimpleHydrogenBond_Dahiyat_Minimization_Creator(double weight , boolean toCalcDerivatives) {
		super(weight);
		this.toCalcDerivatives = toCalcDerivatives;
		parameters = new DahiyatMinimizationParameters();
	}

	public SimpleHydrogenBond_Dahiyat_Minimization_Creator(double weight) {
		this(weight, true);
	}

	public AbstractEnergy createEnergyTerm(Protein protein,
			DistanceMatrix distanceMatrix, CommandList commands) {
		return new SimpleHydrogenBondEnergy(new DahiyatHydrogenBondList(distanceMatrix, protein.atoms(), parameters),
				distanceMatrix,
				weight(),
				toCalcDerivatives);
	}

}
