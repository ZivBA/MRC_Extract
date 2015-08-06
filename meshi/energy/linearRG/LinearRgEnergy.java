package meshi.energy.linearRG;
import java.util.Iterator;

import meshi.energy.CooperativeEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;

public class LinearRgEnergy extends CooperativeEnergyTerm{
    
    public LinearRgEnergy() {}
    
    public LinearRgEnergy(AtomList atomList, 
                    DistanceMatrix dm,
				    double weight) {
	super(toArray(),atomList, dm, null , weight);
	comment = "LinearRG";
    }
 
    public void evaluateAtoms() {
    	double e = evaluate();
    	for (int c=0 ; c<atomList.size() ; c++)
    		atomList.atomAt(c).addEnergy(e/atomList.size());
    }
    
    public double evaluate() {
	if (! on) return 0.0;
 	double energy = 0;
	double Rg = 0;
	double dRg = 0;
	double cmx=0, cmy=0, cmz=0; // center of mass x, y and z
    Atom atom;
	Iterator iter = atomList.iterator(); 
	while((atom = (Atom) iter.next()) != null) { 
	    cmx += atom.x();
	    cmy += atom.y();
	    cmz += atom.z();
	}
	cmx /= atomList.size();
	cmy /= atomList.size();
	cmz /= atomList.size();
	iter = atomList.iterator(); 
	while((atom = (Atom) iter.next()) != null) 
		Rg += (cmx - atom.x())*(cmx - atom.x()) +
		          (cmy - atom.y())*(cmy - atom.y()) +
		          (cmz - atom.z())*(cmz - atom.z()); 
	Rg /= atomList.size();
		
	energy = Math.sqrt(Rg);
	dRg = 1/energy/atomList.size(); // Multiplication by 2 was made to save time
	energy *= weight;
	dRg *= weight;

	iter = atomList.iterator(); 
	while((atom = (Atom) iter.next()) != null) { 
		if (! atom.frozen()) {
		    atom.addToFx(-dRg*(atom.x()-cmx)); // Negating so that it is force
		    atom.addToFy(-dRg*(atom.y()-cmy)); // Negating so that it is force
		    atom.addToFz(-dRg*(atom.z()-cmz)); // Negating so that it is force
		}
    }
    return energy;
    }
}
