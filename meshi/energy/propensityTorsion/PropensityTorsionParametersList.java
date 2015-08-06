package meshi.energy.propensityTorsion;
import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionPair;
import meshi.molecularElements.Residue;
             
public class PropensityTorsionParametersList extends ParametersList {

    public PropensityTorsionParametersList(String[] parametersFileName) {
	    super(parametersFileName,false);  // non-sortable
	}

    public Parameters createParameters(String line) {
    	return new PropensityTorsionParameters(line);	
    }
    
    public Parameters parameters(Object obj) {
	TorsionPair torsionPair = (TorsionPair) obj;
	Torsion torsion1 = torsionPair.torsion1();
	Torsion torsion2 = torsionPair.torsion2();
	Residue residue = torsion1.atom2.residue();
	int resnum1 = torsion1.getTorsionResNum();
	int resnum2 = torsion2.getTorsionResNum();
	int code1 = torsion1.getTorsionCode();
	int code2 = torsion2.getTorsionCode();
	String name1 = torsion1.getTorsionName();
	String name2 = torsion2.getTorsionName();
	PropensityTorsionParameters propensityTorsionParameters;
	if (resnum1 != resnum2) return null; 
	if (code1 < 0) return null; 	
	if (code2 < 0) return null; 	

	for (int cc=0 ; cc<size() ; cc++) {
		propensityTorsionParameters = (PropensityTorsionParameters) elementAt(cc);
    	if ((name1.equals(propensityTorsionParameters.torsion1Name)) &&
	        (name2.equals(propensityTorsionParameters.torsion2Name)) &&
	        (propensityTorsionParameters.mapping[residue.type] > -1)) 
	        return propensityTorsionParameters;
	}

	return null;
    }
}
