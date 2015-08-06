package meshi.energy.twoTorsions;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.ParametersList;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.TorsionPairList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

public class TwoTorsionsCreator extends EnergyCreator  implements KeyWords {
    protected ParametersList parametersList = null;
 
    public TwoTorsionsCreator(double weight) {
	super(weight);
    }
 
    public TwoTorsionsCreator() {
	super(TWO_TORSIONS_ENERGY);
    }
 
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	if (parametersList == null) {
	    // Appanding the path to the filename list of the parameters.
	    String[] strlist = new String[TWO_TORSIONS_PARAMETERS.length];
	    String pathname = parametersDirectory(commands).concat("/");
	    for (int cc=0 ; cc<TWO_TORSIONS_PARAMETERS.length ; cc++)
	        strlist[cc] = pathname.concat(TWO_TORSIONS_PARAMETERS[cc]);
                                 
	    parametersList = new TwoTorsionsParametersList(strlist);
	 }
	TorsionPairList torsionPairList = TorsionPairList.createQuickAndDirtyTorsionPairList(protein,distanceMatrix);
	TorsionPairList relevantTorsionPairList = (TorsionPairList) torsionPairList.filter(new HaveParametersFilter(parametersList),
									 new TorsionPairList()); 
	return new TwoTorsionsEnergy(relevantTorsionPairList, distanceMatrix, (TwoTorsionsParametersList) parametersList, weight(),
				     "2Torsion");
    }

}

