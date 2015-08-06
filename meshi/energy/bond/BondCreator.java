package meshi.energy.bond;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.AtomPairList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;
/**
 * A factory for BondEnergy objects.
 **/
public class BondCreator extends EnergyCreator  implements KeyWords {
    public BondCreator() {
	super(BOND_ENERGY);
    }

    public BondCreator(double weight) {
	super(weight);
    }

    /**
     *<pre> 
     * hides all the hard work needed to generate a BondEnergy object. 
     * a) Extract the bonds from the protein.
     * b) Finds and reads the parameters file.
     **/ 
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
 	AtomPairList bondList = protein.bonds();
	bondList.renumber();
	if (parametersList == null)
	    parametersList = new BondParametersList(parametersDirectory(commands)+
						    "/"+BOND_PARAMETERS);
	return new BondEnergy(bondList, distanceMatrix, (BondParametersList) parametersList, weight());
    }
}
