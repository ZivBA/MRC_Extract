package meshi.energy.cAlphaPlane;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;


public class CAlphaPlaneCreator extends EnergyCreator implements KeyWords {
		
	CAlphaPlaneEnergy planeEnergy;
	
	public CAlphaPlaneCreator() {
		super(CALPHA_HYDROGEN_BONDS_PLANE);    // key for commands file - in the file meshi.util.KeyWords
	}
    
	public CAlphaPlaneCreator(double weight) {
		super(weight);
	}
	
	public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					      CommandList commands) {
        planeEnergy = new CAlphaPlaneEnergy(protein, distanceMatrix,
                               null,weight());
		return planeEnergy;
	}
		

	public final CAlphaPlaneEnergy getPlaneEnergy() {
		return planeEnergy;
	}
}
