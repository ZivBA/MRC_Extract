package meshi.energy.bond;
import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.AtomPair;

public class BondEnergyElement extends EnergyElement {
    protected Atom atom1, atom2;
    protected int number1, number2;
    protected Distance distance;
    protected double target, force, force2;
    protected boolean frozen;
    double weight;
    public BondEnergyElement() {}
    public BondEnergyElement(AtomPair atomPair, Parameters parameters, 
			     DistanceMatrix distanceMatrix, double weight) {
	this.weight = weight;
	atom1 = atomPair.atom1();
	atom2 = atomPair.atom2();
	setAtoms();
	int atom1Number = atomPair.atom1Number();
	int atom2Number = atomPair.atom2Number();
	distance = distanceMatrix.distance(atom1Number, atom2Number);			   
	target = ((BondParameters) parameters).target;
	force = ((BondParameters) parameters).force*weight;
	force2 = ((BondParameters) parameters).force2*weight;	    
	updateFrozen();
    }
    protected void setAtoms(){
	atoms = new AtomList();
	atoms.add(atom1);
	atoms.add(atom2);
    }
	
	

    public double evaluate() {
	double d,d2,d3;
	double deDd;
	double deDx;
	double deDy;
	double deDz;
	double energy;
	double ALPHA = 0.0001;
	double d2PlusAlpha;
	double d2PlusAlpha2;
	double dis = distance.distance();
	if (frozen()) return 0.0;
	d = dis - target;
	d2 = d*d;
	energy = d2 * force;
	deDd =  d * force2;

	if (dis < 0.2) {
		energy += (dis-0.2)*(dis-0.2)*100000;
		deDd += 200000*(dis-0.2);
	}
//	d3 = d2*d;
//	d2PlusAlpha = d2+ALPHA;
//	d2PlusAlpha2 = d2PlusAlpha*d2PlusAlpha;
//	energy += force*d2/d2PlusAlpha;
//	double deDdTemp = 2*force*(d/d2PlusAlpha - d3/d2PlusAlpha2);
//	deDd += deDdTemp;

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
	//System.out.println(energy +  " " + atom1.residueNumber() + " " + atom1.name() + " " + atom2.residueNumber() + " " + atom2.name());
	return energy;
    }

    public String toString() {
	return ("BondEnergyElement target = "+dFormatSrt.f(target)+" force = "+dFormatSrt.f(force)+" distance = "+
		dFormatSrt.f(distance.distance())+" energy = "+dFormatSrt.f(evaluate())+"\n"+
		atom1.verbose(1)+"\n"+atom2.verbose(1));
    }
}
