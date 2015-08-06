package meshi.energy.inflate;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

public class InflateCreator extends EnergyCreator  implements KeyWords {
	private double rmsTarget = 0.0;

    public InflateCreator(double weight, double rmsTarget) {
    	super(weight);
    	this.rmsTarget = rmsTarget;
    }

    public InflateCreator() {
	super(0.0);
    }
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	return new Inflate(distanceMatrix, weight(), rmsTarget);
    }
}
