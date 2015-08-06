package meshi.energy.CAsolvate;
import meshi.energy.Parameters;
import meshi.energy.ParametersList;

/**
 * The CAsolvate class does not need any paramters from files. It get them all in the creator.
 **/

public class CAsolvateParametersList extends ParametersList {
	
    public CAsolvateParametersList(String[] parameterFiles) {
	super();
	}
    
    public Parameters createParameters(String s) {
    	throw new RuntimeException("This method should not be called");
    }
    
    public Parameters parameters(Object obj) {
    	    	throw new RuntimeException("This method should not be called");
    }
 
}
