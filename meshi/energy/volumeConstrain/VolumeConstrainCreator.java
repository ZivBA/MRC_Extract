package meshi.energy.volumeConstrain;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

public class VolumeConstrainCreator extends EnergyCreator implements KeyWords {
    public VolumeConstrainCreator() {
	super(VOLUME_CONSTRAIN);
    }
    
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	throw new RuntimeException("not implemented");
    }
}
