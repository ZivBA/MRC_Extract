package meshi.energy.twoTorsions;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.TorsionPairList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

public class FlatRamachCreator extends EnergyCreator  implements KeyWords {

    public FlatRamachCreator(double weight) {
	super(weight);
    }

    public FlatRamachCreator() {
	super(FLAT_RAMACH_ENERGY);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	if (parametersList== null) {
	    // Appanding the path to the filename list of the parameters.
	    String[] strlist = new String[FLAT_RAMACH_PARAMETERS.length];
	    String pathname = parametersDirectory(commands).concat("/");
	    for (int cc=0 ; cc<FLAT_RAMACH_PARAMETERS.length ; cc++)
	        strlist[cc] = pathname.concat(FLAT_RAMACH_PARAMETERS[cc]);
                                 
	    parametersList = new TwoTorsionsParametersList(strlist);
	 }

	TorsionPairList torsionPairList = TorsionPairList.createQuickAndDirtyTorsionPairList(protein,distanceMatrix);
	TorsionPairList relevantTorsionPairList = (TorsionPairList) torsionPairList.filter(new HaveParametersFilter(parametersList),
									 new TorsionPairList()); 

	return new TwoTorsionsEnergy(relevantTorsionPairList, distanceMatrix, (TwoTorsionsParametersList) parametersList, weight(),
				     "FlatRamachandran");
    }
}
