package meshi.energy.cAlphaHydrogenBond;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;


public class CAlphaHydrogenBondsCreator extends EnergyCreator implements KeyWords {
		
	CAlphaHydrogenBondsEnergy hydrogenBondsEnergy;
	
	public CAlphaHydrogenBondsCreator() {
		super(CALPHA_HYDROGEN_BONDS);    // key for commands file - in the file meshi.util.KeyWords
	}
    
	public CAlphaHydrogenBondsCreator(double weight) {
		super(weight);
	}
	
	public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					      CommandList commands) {	        
        hydrogenBondsEnergy = new CAlphaHydrogenBondsEnergy(protein, distanceMatrix,
                               null,weight());
		return hydrogenBondsEnergy;
	}
		

	public final CAlphaHydrogenBondsEnergy getHydrogenBondsEnergy() {
		return hydrogenBondsEnergy;
	}
}
