package meshi.energy.alphaTorsion;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.AngleList;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.TorsionList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

public class AlphaTorsionCreator extends EnergyCreator  implements KeyWords {
    public AlphaTorsionCreator(double weight) {
	super(weight);
    }
    public AlphaTorsionCreator() {
	super(ALPHA_TORSION_ENERGY);
    }
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
					  	
	if (parametersList== null) {                                 
	    parametersList = new AlphaTorsionParametersList(parametersDirectory(commands)+
						     "/"+ALPHA_TORSION_PARAMETERS);
	}
	AngleList angleList = AngleList.getCaAngles(protein,distanceMatrix).namedFilter();
	TorsionList torsionList = new TorsionList(angleList, distanceMatrix);
	TorsionList relevantTorsionList = (TorsionList)torsionList.filter(new HaveParametersFilter(parametersList),
									 new TorsionList()); 

	return new AlphaTorsionEnergy(relevantTorsionList, distanceMatrix, (AlphaTorsionParametersList) parametersList, weight());
    }
}
