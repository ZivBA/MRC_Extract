package meshi.energy.ROT1solvation.parameters;

import meshi.geometry.Distance;

public abstract class AbstractCBParameters extends AbstractROT1Parameters {

	public AbstractCBParameters(String parameterFileName) {
		super(parameterFileName);
	}

// the code for CA and CB representaion:
//	public boolean filterDisForRelevance(Distance dis) {
//		if (!((dis.atom1().name().equals("CB") || dis.atom1().name().equals("CA")) && 
//				(dis.atom2().name().equals("CB") || dis.atom2().name().equals("CA"))))
//			return false;
//		if ((dis.atom1().residueNumber()>(dis.atom2().residueNumber()-2)) && 
//				(dis.atom1().residueNumber()<(dis.atom2().residueNumber()+2)))
//			return false;
//		return true;
//	}

	public boolean filterDisForRelevance(Distance dis) {
		if (!dis.atom1().name().equals("CB") || !dis.atom2().name().equals("CB"))
			return false;
		if ((dis.atom1().residueNumber()>(dis.atom2().residueNumber()-2)) && 
				(dis.atom1().residueNumber()<(dis.atom2().residueNumber()+2)))
			return false;
		return true;
	}

}
