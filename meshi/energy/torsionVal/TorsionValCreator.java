package meshi.energy.torsionVal;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.TorsionList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

           
public class TorsionValCreator extends EnergyCreator  implements KeyWords {

    public TorsionValCreator(double weight) {
	super(weight);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	if (parametersList== null)
	    parametersList = new TorsionValParametersList();
	TorsionList torsionList = (TorsionList) (TorsionList.createTorsionList(protein, distanceMatrix)).filter(
		new TorsionList.FilterFamouseTorsions(),new TorsionList());
	return new TorsionValEnergy(distanceMatrix, torsionList, (TorsionValParametersList) parametersList, weight());
    }
}
