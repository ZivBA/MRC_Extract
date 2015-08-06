package meshi.energy.disConst;

import meshi.energy.EnergyElement;
import meshi.geometry.Distance;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;


public class DisConstEnergyElement extends EnergyElement {
    protected double weight;
    protected Distance distance;
    protected double targetDis;
    protected Atom atom1;
    protected Atom atom2;
    
    public DisConstEnergyElement() {}
    
    public DisConstEnergyElement(double weight, Distance dis) {
        this.weight = weight;
        distance = dis;
        atom1 = distance.atom1();
        atom2 = distance.atom2();
        targetDis = dis.distance();
        setAtoms();
        updateFrozen();
    }

    protected void setAtoms(){
    	atoms = new AtomList();
    	atoms.add(atom1);
    	atoms.add(atom2);
    }

    public void setTargetDis(double newTarget) {
    	targetDis = newTarget;
    }
    
    public double evaluate() {
	double d;
	double deDd;
	double deDx;
	double deDy;
	double deDz;
	double energy;

    distance.update();
    d = (distance.distance()-targetDis) ;
    if (frozen() || d==0.0) return 0.0;

    energy = d * d * weight;
	deDd =  2 * d * weight;

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

    public void update() {}
                                          
    public String toString() {
    	return atom1.toString()+"\n"+atom2.toString()+"\n";
    }
}
