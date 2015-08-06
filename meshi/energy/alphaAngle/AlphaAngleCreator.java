package meshi.energy.alphaAngle;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.AngleList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

public class AlphaAngleCreator extends EnergyCreator  implements KeyWords {
    public AlphaAngleCreator(double weight) {
	super(weight);
    }
    public AlphaAngleCreator() {
	super(ALPHA_ANGLE_ENERGY);
    }
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	AngleList angleList = AngleList.getCaAngles(protein,distanceMatrix).namedFilter();
	if (parametersList== null) {                                 
	    parametersList = new AlphaAngleParametersList(parametersDirectory(commands)+
						     "/"+ALPHA_ANGLE_PARAMETERS);
	 }
	return new AlphaAngleEnergy(angleList, distanceMatrix, (AlphaAngleParametersList) parametersList, weight());
    }
}
