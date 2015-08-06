package meshi.energy.outOfPlane;
import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.geometry.Torsion;

public class OutOfPlaneParametersList extends ParametersList {
    public OutOfPlaneParametersList(String parametersFileName) {
	super(parametersFileName,true);
    }

    public Parameters createParameters(String line) {
	return new OutOfPlaneParameters(line);
    }

    public Parameters parameters(Object obj) {
	Torsion torsion = (Torsion) obj;
	Parameters key = new OutOfPlaneParameters(torsion.atom1.type, torsion.atom2.type,
						  torsion.atom3.type, torsion.atom4.type);
	return getParameters(key);
    }
}
