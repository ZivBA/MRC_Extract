
package meshi.energy.hydrogenBondsPlane;

import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.energy.ParametersList;


public class HydrogenBondsPlaneParametersList  extends ParametersList {

	public HydrogenBondsPlaneParametersList(String parametersFileName){
		super(parametersFileName, false);                                           
	}

        public Parameters createParameters(String line) {
		return new HydrogenBondsPlaneParameters (new StringTokenizer(line));
	}

   public Parameters parameters(Object obj) {
    	return null;
    }

}
