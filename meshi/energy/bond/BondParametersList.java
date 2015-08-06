package meshi.energy.bond;
import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.molecularElements.AtomPair;
/**
 * ParametersList for BondEnergy. The list is sortable.
 **/
public class BondParametersList extends ParametersList {
        public BondParametersList(String parametersFileName) {
	    super(parametersFileName,true);
	}

     public Parameters createParameters(String line) {
	return new BondParameters(line);
    }
    
    /**
     * Returns BondParameters object (which specify target-distance and force-constant) 
     * for a pair of atoms).
     **/ 
    public Parameters parameters(Object obj) {
	AtomPair pair = (AtomPair) obj;
	Parameters key = new BondParameters(pair.largeType(), pair.smallType());
	return getParameters(key);
    }
}
