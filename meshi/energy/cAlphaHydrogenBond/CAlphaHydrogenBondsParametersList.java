package meshi.energy.cAlphaHydrogenBond;

import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.energy.ParametersList;

public class CAlphaHydrogenBondsParametersList extends ParametersList {

	public CAlphaHydrogenBondsParametersList(String parametersFileName){
		super(parametersFileName, false);
	}

	public Parameters createParameters(String line) {
			return new CAlphaHydrogenBondsParameters(new StringTokenizer(line));
	}
        
    public Parameters parameters(Object obj) {
	return null;
    }
}
