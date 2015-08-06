package meshi.energy.LennardJones;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

public class LennardJonesCreator extends EnergyCreator  implements KeyWords {


	private int type=0; // allows us to run one of the EV scenario 
	private double frac=1.0; // allows us to multiply the atomic radii by a certain fraction 
	
	

    public LennardJonesCreator(double weight , int type , double frac) {
  	super(weight);
  	this.type = type;
  	this.frac = frac;
    }
 
    public LennardJonesCreator(double weight , int type) {
  	super(weight);
  	this.type = type;
    }

    public LennardJonesCreator(double weight) {
  	super(weight);
    }

    public LennardJonesCreator() {
  	super(EXCLUDED_VOL);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	if (parametersList== null)
	    parametersList = new LennardJonesParametersList(parametersDirectory(commands)+
							    "/"+EXCLUDED_VOL_PARAMETERS , frac);
	return new LennardJones(distanceMatrix, ( LennardJonesParametersList) parametersList, type,  weight());
    }
}
