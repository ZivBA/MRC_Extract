package meshi.geometry;
import meshi.molecularElements.Atom;

/**
 *
 *--------------------------- TorsionForAnyDistance -----------------------------------
 *TorsionForAnyDistance is equal to Torsion, but it calculates only one predefined type of torsions.
 *Besides, it checks  if some distance does not exist in DistanceMatrix than this distance will be calculated here
 **/
	
public class DisposableTorsion extends Torsion{

    public DisposableTorsion(Atom atom1, Atom atom2, Atom atom3, Atom atom4, DistanceMatrix distanceMatrix) {
     super(atom1,atom2,atom3,atom4,getDistance(atom1,atom2,distanceMatrix),
     		getDistance(atom3,atom2,distanceMatrix),getDistance(atom3,atom4,distanceMatrix));   
    }

    public void update(int numberOfUpdates){
		throw new RuntimeException("Weird using of DisposableTorsion.update\n"+
				"DisposableTorsion is not updatable.");
	}
    
	private static Distance getDistance(Atom atom1, Atom atom2, DistanceMatrix distanceMatrix){
        Distance out = distanceMatrix.distance(atom1,atom2);
        if (out.distance() == Distance.INFINITE_DISTANCE)
                out = new Distance (atom1, atom2);
		return out;
	}    
}
