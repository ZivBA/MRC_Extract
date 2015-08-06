package meshi.energy.simpleHydrogenBond;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.hydrogenBond.Dahiyat.DahiyatHydrogenBondListBBonlyInLoop;
import meshi.geometry.hydrogenBond.Dahiyat.DahiyatLowAccuracyParameters;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;

public class SimpleHydrogenBond_Dahiyat_LowAccuracy_BBonly_InLoop_Creator extends EnergyCreator {

	private int loopStart=-1;
	private int loopEnd=-1;
	private boolean toCalcDerivatives = true;
	private DahiyatLowAccuracyParameters parameters = null;

	public SimpleHydrogenBond_Dahiyat_LowAccuracy_BBonly_InLoop_Creator(double weight , boolean toCalcDerivatives, 
			int loopStart, int loopEnd) {
		super(weight);
		this.toCalcDerivatives = toCalcDerivatives;
		parameters = new DahiyatLowAccuracyParameters();
		this.loopStart = loopStart;
		this.loopEnd = loopEnd;
	}

	public SimpleHydrogenBond_Dahiyat_LowAccuracy_BBonly_InLoop_Creator(double weight,int loopStart, int loopEnd) {
		this(weight, true, loopStart, loopEnd);
	}

	public AbstractEnergy createEnergyTerm(Protein protein,
			DistanceMatrix distanceMatrix, CommandList commands) {
		return new SimpleHydrogenBondEnergy(new DahiyatHydrogenBondListBBonlyInLoop(distanceMatrix, protein.atoms(), parameters, loopStart, loopEnd),
				distanceMatrix,
				weight(),
				toCalcDerivatives);
	}

}
