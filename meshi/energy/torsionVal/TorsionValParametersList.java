package meshi.energy.torsionVal;
import meshi.energy.Parameters;
import meshi.energy.ParametersList;

public class TorsionValParametersList extends ParametersList {
    public TorsionValParametersList() {
	super();
    }

    public Parameters createParameters(String line) {
	return new TorsionValParameters();
    }

    public Parameters parameters(Object obj) {
	return new TorsionValParameters();
    }
}
