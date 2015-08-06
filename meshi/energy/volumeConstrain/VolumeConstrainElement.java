package meshi.energy.volumeConstrain;
import meshi.energy.EnergyElement;
import meshi.geometry.Distance;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;

public class VolumeConstrainElement extends EnergyElement {
    Distance distance;
    double radius, weight;
    Atom atom;
    public VolumeConstrainElement() {}
    public VolumeConstrainElement(Atom center, Atom atom, double radius, double weight) {
	distance = new Distance(center,atom);
	this.radius = radius;
	this.weight = weight;
	this.atom = atom;
	setAtoms();
    }

    public void setAtoms() {
	atoms = new AtomList();
	atoms.add(atom);
    }
	
    public double evaluateAtoms() {return 0.0;}

    public double evaluate() {
	double energy = 0;
	double d,force;
	
	distance.update();
	d = distance.distance() - radius;
	if (d > 0) { 
	    // that is the atom is out of the sphere
	    energy = weight * d*d;
	// The force is the negative of the derivative.
	    // Note that the distance derivatives are reported for the first point,
	    // and we use the second (the atom).
	    force = 2.0 * weight * d; 
	    atom.addToFx(force * distance.dDistanceDx());
	    atom.addToFy(force * distance.dDistanceDy());
	    atom.addToFz(force * distance.dDistanceDz());
	}
	return energy;
    }
	
}
    
