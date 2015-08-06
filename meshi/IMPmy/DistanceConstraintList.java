package meshi.IMPmy;

import java.util.Vector;

public class DistanceConstraintList extends Vector<DistanceConstraint> {
	
	public String toString() {
		String out = "";
		for (DistanceConstraint disConst : this) {
//			if ((disConst.constraintType==DistanceConstraintType.CROSS_LINK) || (disConst.constraintType==DistanceConstraintType.CONNECTIVITY)) {
			if (disConst.evaluate(false)>0.2) {
				out += (disConst + "\n");
			}
		}
		return out;
	}

	public void zeroConnectivity() {
		for (DistanceConstraint dis : this) {
			if (dis.constraintType == DistanceConstraintType.CONNECTIVITY) {
				dis.setWeight(0.0);
			}			
		}		
	}
	
	public void upRigidity() {
		for (DistanceConstraint dis : this) {
			if (dis.constraintType == DistanceConstraintType.RIGID_BODY) {
				dis.setWeight(1.0);
			}			
		}		
	}
	
}
