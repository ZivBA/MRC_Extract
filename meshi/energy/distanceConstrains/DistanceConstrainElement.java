package meshi.energy.distanceConstrains;
import meshi.energy.EnergyElement;
import meshi.geometry.Distance;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
//------------------------------------------------------------------------------------------------
public class DistanceConstrainElement extends EnergyElement {
    protected Atom atom1, atom2;
    protected int number1, number2;
    protected Distance distance;
    protected double target, tolerance, force;
    protected double weight;

    public DistanceConstrainElement(Atom atom1, Atom atom2, DistanceConstrainParameters parameters, 
			     double weight) {
	this.weight = weight;
	this.atom1 = atom1;
	this.atom2 = atom2;
	setAtoms();
	number1 = atom1.number();
	number2 = atom2.number();
	distance = new Distance(atom1, atom2);
	target = parameters.target;
	force =  parameters.force*weight;
	tolerance =  parameters.tolerance;
	updateFrozen();
    }
	
    protected void setAtoms(){
	atoms = new AtomList();
	atoms.add(atom1);
	atoms.add(atom2);
    }

    public double evaluate() {
	double dis;
	double deDd = 0;
	double deDx;
	double deDy;
	double deDz;
	double energy = 0;

	if (frozen()) return 0.0;
	distance.update();
	dis = distance.distance();
	if (dis >= (target + tolerance)) {
		energy = force*(dis - (target + tolerance))*(dis - (target + tolerance));
		deDd = 2*force*(dis - (target + tolerance));
	}
	else if (dis <= (target - tolerance)) {
		energy = force*(dis - (target - tolerance))*(dis - (target - tolerance));
		deDd = 2*force*(dis - (target - tolerance));
	}
	else {
		energy = 0.0;
		deDd = 0.0;
	}

	deDx = deDd*distance.dDistanceDx();
	deDy = deDd*distance.dDistanceDy();
	deDz = deDd*distance.dDistanceDz();
	if (! atom1.frozen()) {
	    atom1.addToFx(-1*deDx); // force = -derivative   
	    atom1.addToFy(-1*deDy);
	    atom1.addToFz(-1*deDz);
	}
	if (! atom2.frozen()) {
	    atom2.addToFx(deDx);
	    atom2.addToFy(deDy);
	    atom2.addToFz(deDz);
	}
	return energy;
    }
    
    public String toString() {
	return ("DistanceConstrain target = "+dFormatSrt.f(target)+" force = "+dFormatSrt.f(force)+" distance = "+
		dFormatSrt.f(distance.distance())+" energy = "+dFormatSrt.f(evaluate())+" "+
		"tolerance "+dFormatSrt.f(tolerance)+"\n"+
		atom1.verbose(1)+"\n"+atom2.verbose(1));
    }
	
}
