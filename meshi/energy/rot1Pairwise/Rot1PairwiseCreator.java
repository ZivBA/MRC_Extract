package meshi.energy.rot1Pairwise;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

public class Rot1PairwiseCreator extends EnergyCreator  implements KeyWords {


	private int type=0; // allows us to run one of the EV scenario 
	private double frac=1.0; // allows us to multiply the atomic radii by a certain fraction 
	

    public Rot1PairwiseCreator(double weight , int type , double frac) {
  	super(weight);
  	this.type = type;
  	this.frac = frac;
    }
 
    public Rot1PairwiseCreator(double weight , int type) {
  	super(weight);
  	this.type = type;
    }

    public Rot1PairwiseCreator(double weight) {
  	super(weight);
    }

    public Rot1PairwiseCreator() {
  	super(EXCLUDED_VOL);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	if (parametersList== null) 
	    parametersList = new Rot1PairwiseParametersList(parametersDirectory(commands)+
			    "/"+ROT1_PAIRWISE_PARAMETERS);
	return new Rot1Pairwise(distanceMatrix, ( Rot1PairwiseParametersList) parametersList, type,  weight());
    }
}
