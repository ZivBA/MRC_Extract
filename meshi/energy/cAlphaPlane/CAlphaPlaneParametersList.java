package meshi.energy.cAlphaPlane;

import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.energy.ParametersList;

public class CAlphaPlaneParametersList extends ParametersList {

	public CAlphaPlaneParametersList(String parametersFileName){
		super(parametersFileName, false);
	}

        public Parameters createParameters(String line) {
		return new CAlphaPlaneParameters(new StringTokenizer(line));
	}
        
public Parameters parameters(Object obj) {
	return null;
    }

}
