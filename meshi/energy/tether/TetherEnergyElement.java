package meshi.energy.tether;

import meshi.energy.EnergyElement;
import meshi.geometry.Distance;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;


public class TetherEnergyElement extends EnergyElement {
    protected Atom atom;   
    protected double weight;
    protected boolean frozen;
    protected Distance distance;
    protected Atom  evaluatedAtom;
    
    public TetherEnergyElement() {}
    
    public TetherEnergyElement(double weight, Atom inputAtom, Atom pegAtom) {
        atom = inputAtom;
        
        this.weight = weight;
//        if (atom.isBackbone)
//        	this.weight = 1*weight;
//        else
//        	this.weight = 0.15*weight;
        
        
        setAtoms();
        updateFrozen();

        if (pegAtom==null)
        	evaluatedAtom =  new Atom(atom.x(),atom.y(),atom.z(),"eval","dummy",1,2);
        else
        	evaluatedAtom = new Atom(pegAtom.x(),pegAtom.y(),pegAtom.z(),"eval","dummy",1,2);
        distance = new Distance(atom,evaluatedAtom);
    }

    protected void setNewPegCoordinates(double newX,double newY,double newZ) {
    	evaluatedAtom.setXYZ(newX, newY, newZ);
    }
    
    
    protected void setAtoms(){
    	atoms = new AtomList();;
    	atoms.add(atom);
    }

    public double evaluate() {
   	if (frozen() || (atom.reliability()==0.0)) return 0.0;    	
	double d;
	double deDd;
	double deDx;
	double deDy;
	double deDz;
	double energy;

    distance.update();
    d = distance.distance() ;
    if (d==0.0) return 0.0;

    energy = d * d * weight * atom.reliability();
	deDd =  2 * d * weight * atom.reliability();



	deDx = deDd*distance.dDistanceDx();
	deDy = deDd*distance.dDistanceDy();
	deDz = deDd*distance.dDistanceDz();
	if (! atom.frozen()) {
	    atom.addToFx(-1*deDx); // force = -derivative
	    atom.addToFy(-1*deDy);
	    atom.addToFz(-1*deDz);
	}
	
	return energy;
    }

    public void update() {}
                                          
    public String toString() {
    	return atom.toString();
    }
}
