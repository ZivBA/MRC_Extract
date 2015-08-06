package meshi.energy.angle;
import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.geometry.Angle;

public class AngleParametersList extends ParametersList {
        public AngleParametersList(String parametersFileName) {
	    super(parametersFileName,true);
	}

     public Parameters createParameters(String line) {
	return new AngleParameters(line);
    }

    public Parameters parameters(Object obj) {
	Angle angle = (Angle) obj;
	Parameters key = new AngleParameters(angle.atom1.type, angle.atom2.type, angle.atom3.type);
	return getParameters(key);
    }
}
