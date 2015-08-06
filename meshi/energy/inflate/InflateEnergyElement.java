package meshi.energy.inflate;
import meshi.energy.EnergyElement;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.MeshiProgram;
//-------------------------------------------------------------------------------------------------------
//                                     InflateElement
//-------------------------------------------------------------------------------------------------------
public  class InflateEnergyElement extends EnergyElement {
    protected Atom atom1, atom2;
    //    protected AtomPair atomPair;
    protected Atom atom1Copy, atom2Copy;
    protected double dEdD;
    protected double dEdX;
    protected double dEdY;
    protected double dEdZ;
    protected double energy;
    protected double weight;
    protected Distance distance;
    protected double target;
    public static final double ALPHA = 0.1;
    private double rMax2;
    
    public  InflateEnergyElement(Distance atoms,  DistanceMatrix distanceMatrix, double weight) {
	atom1 = atoms.atom1();
	atom2 = atoms.atom2();
        //atomPair = atoms;
	setAtoms();
	updateFrozen();
        rMax2 = DistanceMatrix.rMax2();
	distance = atom1.distance(atom2);
	this.weight = MeshiProgram.randomNumberGenerator().nextDouble()* weight; 
        target = distanceMatrix.rMax();
    }

    protected void setAtoms(){
	atoms = new AtomList();
	atoms.add(atom1);
	atoms.add(atom2);
    }
	    	
	
    public double evaluate() {
	double dis;
	double d, d2, twoD, invDplus, invDplus2;
	double deDd = 0;
	double deDx;
	double deDy;
	double deDz;
	double energy = 0;

	if (frozen()) return 0;
	distance.update();
	if (distance.distance2()>rMax2) return 0;
	dis = distance.distance();
	d = dis - target;
	energy = -1 * weight*d*d/(d-ALPHA);
	deDd   = -1 * weight*(1 - ALPHA*ALPHA/((d-ALPHA)*(d-ALPHA)));

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
}
