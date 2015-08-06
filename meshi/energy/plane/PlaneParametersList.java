
package meshi.energy.plane;
import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.geometry.Torsion;
             
public class PlaneParametersList extends ParametersList {
        public PlaneParametersList(String parametersFileName) {
	    super(parametersFileName,true);
	}

     public Parameters createParameters(String line) {
	return new PlaneParameters(line);
    }
    public Parameters parameters(Object obj) {
	Torsion torsion = (Torsion) obj;
	Parameters key = new PlaneParameters(torsion.atom1.type, torsion.atom2.type,
					     torsion.atom3.type, torsion.atom4.type);
	return getParameters(key);
    }
}
