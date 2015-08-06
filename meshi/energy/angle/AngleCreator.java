package meshi.energy.angle;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.AngleList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.AtomPairList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

public class AngleCreator extends EnergyCreator  implements KeyWords {
    public AngleCreator() {
	super(ANGLE_ENERGY);
    }

    public AngleCreator(double weight) {
	super(weight);
    }
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
 	AtomPairList bondList = protein.bonds();
	bondList.renumber();
	AngleList angleList = new AngleList(bondList, distanceMatrix);
	if (parametersList== null)
	    parametersList = new AngleParametersList(parametersDirectory(commands)+
						     "/"+ANGLE_PARAMETERS);
	return new AngleEnergy(angleList, distanceMatrix, (AngleParametersList) parametersList, weight());
    }
}
