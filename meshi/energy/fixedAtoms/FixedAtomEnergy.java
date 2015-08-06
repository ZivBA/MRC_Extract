package meshi.energy.fixedAtoms;
import meshi.energy.CooperativeEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;


/**
 **/

public class FixedAtomEnergy extends CooperativeEnergyTerm{

    private double wid,wid22,invWid22,cutoff,uprange;
    private double[][] aux;
    private int[] lut; 
    private int maxAtomNum;
    private Atom atom;
    private double[] parameterX, parameterY, parameterZ; 

    public FixedAtomEnergy() {}
    
    public FixedAtomEnergy(AtomList atomList, 
                    DistanceMatrix dm,
				    double weight) {
	super(toArray(),atomList, dm, null, weight);
	comment = "FixedAtoms";
	parameterX = new double[atomList.size()];
	parameterY = new double[atomList.size()];
	parameterZ = new double[atomList.size()];
	for (int c=0; c<atomList.size() ; c++) {
	   atom = atomList.atomAt(c);
	   parameterX[c] = atom.x();	
	   parameterY[c] = atom.y();	
	   parameterZ[c] = atom.z();	
	}
    }
    
    public void evaluateAtoms() {
	      evaluate(true);
    }
    
    public double evaluate() {
    	double e = evaluate(false);
    	return e;
    }

    public double evaluate(boolean updateAtoms) {
	double energy = 0;
	int c;
       
	//Reseting the aux matrix
	for (c=0 ; c<atomList.size() ; c++) {
	    atom = atomList.atomAt(c);
	    energy += weight*((atom.x()-parameterX[c])*(atom.x()-parameterX[c]) +
   	                     (atom.y()-parameterY[c])*(atom.y()-parameterY[c]) +
   	                     (atom.z()-parameterZ[c])*(atom.z()-parameterZ[c]));
   	    if (updateAtoms)
   	        atom.addEnergy(weight*((atom.x()-parameterX[c])*(atom.x()-parameterX[c]) +
   	                     (atom.y()-parameterY[c])*(atom.y()-parameterY[c]) +
   	                     (atom.z()-parameterZ[c])*(atom.z()-parameterZ[c])));
		if (! atom.frozen()) {
		    atom.addToFx(-weight*2*(atom.x()-parameterX[c]));
		    atom.addToFy(-weight*2*(atom.y()-parameterY[c]));
		    atom.addToFz(-weight*2*(atom.z()-parameterZ[c]));
		}
	}        
    return energy;
    }
    
}
	
	
