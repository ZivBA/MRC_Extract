package meshi.energy.plane;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.TorsionList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

public class PlaneCreator extends EnergyCreator  implements KeyWords {
    public PlaneCreator(double weight) {
	super(weight);
    }
    public PlaneCreator() {
	super(PLANE_ENERGY);
    }
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	if (parametersList== null)
	    parametersList = new PlaneParametersList(parametersDirectory(commands)+
						     "/"+PLANE_PARAMETERS);
	
	TorsionList torsionList = TorsionList.createQuickAndDirtyTorsionList(protein, distanceMatrix);
	TorsionList relevantTorsionList = (TorsionList)torsionList.filter(new HaveParametersFilter(parametersList),
									 new TorsionList());
	return new PlaneEnergy(relevantTorsionList, distanceMatrix, (PlaneParametersList) parametersList, weight());
    }
}
	    
