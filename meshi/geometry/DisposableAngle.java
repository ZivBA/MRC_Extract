package meshi.geometry;
import meshi.molecularElements.Atom;
import meshi.util.UpdateableException;
/**
 *
 *--------------------------- Angle -----------------------------------
 **/
	
public class DisposableAngle  extends Angle {
    public DisposableAngle(Atom atom1, Atom atom2, Atom atom3, DistanceMatrix distanceMatrix) {
		super(atom1, atom2, atom3, 
			   getDistance(atom1,atom2,distanceMatrix),
			   getDistance(atom3,atom2,distanceMatrix));
    }

    public void update(int numberOfUpdates) throws UpdateableException{
		throw new RuntimeException("Weird using of DisposableAngle.update\n"+
				"DisposableAngle is not updatable.");
	}


	private static Distance getDistance(Atom atom1, Atom atom2, DistanceMatrix distanceMatrix){
        Distance out = distanceMatrix.distance(atom1,atom2);
        if (out.distance() == Distance.INFINITE_DISTANCE)
                out = new Distance (atom1, atom2);
		return out;
	}
}
