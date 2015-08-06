package meshi.energy.propensityTorsion;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.TorsionPairList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

public class PropensityTorsionCreator extends EnergyCreator  implements KeyWords {
    public PropensityTorsionCreator(double weight) {
	super(weight);
    }
    public PropensityTorsionCreator() {
	super(PROPENSITY_TORSION_ENERGY);
    }
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	if (parametersList== null) {
	    // Appanding the path to the filename list of the parameters.
	    String[] strlist = new String[PROPENSITY_TORSION_PARAMETERS.length];
	    String pathname = parametersDirectory(commands).concat("/");
	    for (int cc=0 ; cc<PROPENSITY_TORSION_PARAMETERS.length ; cc++)
	        strlist[cc] = pathname.concat(PROPENSITY_TORSION_PARAMETERS[cc]);
                                 
	    parametersList = new PropensityTorsionParametersList(strlist);
	 }

	TorsionPairList torsionPairList = TorsionPairList.createTorsionPairList(protein,distanceMatrix);
	TorsionPairList relevantTorsionPairList = (TorsionPairList) torsionPairList.filter(new HaveParametersFilter(parametersList),
									 new TorsionPairList()); 

	return new PropensityTorsionEnergy(relevantTorsionPairList, distanceMatrix, (PropensityTorsionParametersList) parametersList, weight());

    }
}
