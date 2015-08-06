package meshi.applications.loopBuilding;

import java.util.Vector;

public class BasicLoopResultVector extends Vector<BasicLoopResult> {

private static final long serialVersionUID = 1L;

public boolean hasNoSimilar(double[][] refLoopCoors, double rmsCriterion) {
	if (rmsCriterion<-0.1)
		return true;
	for (int c=0 ; c<size() ; c++)
		if (get(c).calcRMS(refLoopCoors)<rmsCriterion)
			return false;
	return true;
}

}
