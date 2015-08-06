package meshi.energy.outOfPlane;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.TorsionList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

           
public class OutOfPlaneCreator extends EnergyCreator  implements KeyWords {
    public OutOfPlaneCreator(double weight) {
	super(weight);
    }
    public OutOfPlaneCreator() {
	super(OUT_OFPLANE_ENERGY);
    }
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	if (parametersList== null)
	    parametersList = new OutOfPlaneParametersList(parametersDirectory(commands)+
							  "/"+OUT_OF_PLANE_PARAMETERS);
	TorsionList torsionList = TorsionList.createQuickAndDirtyTorsionList(protein, distanceMatrix);
	    TorsionList OOPlist = (TorsionList) torsionList.filter(new TorsionList.FilterOOP(), 
                                      new TorsionList());
	TorsionList relevantTorsionList = (TorsionList) OOPlist.filter(new HaveParametersFilter(parametersList),
									  new TorsionList());
	return new OutOfPlaneEnergy(distanceMatrix, relevantTorsionList, (OutOfPlaneParametersList) parametersList, weight());
    }
}
