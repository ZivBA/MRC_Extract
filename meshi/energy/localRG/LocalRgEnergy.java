package meshi.energy.localRG;
import meshi.energy.CooperativeEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;

public class LocalRgEnergy extends CooperativeEnergyTerm{
    
   	private int begins=-1,ends=-1;
   	private double[] distances;
   	private Atom[] atoms;
    
    public LocalRgEnergy() {}
    
    public LocalRgEnergy(AtomList atomList, 
                    DistanceMatrix dm,
                    int begins,
                    int ends,
				    double weight) {
	super(toArray(dm),atomList, dm, null , weight);
	this.begins = begins;
	this.ends = ends;
	distances = new double[ends-begins+1];
	atoms = new Atom[ends-begins+1];
	int cc=0;	
	for (int c=0 ; c<atomList.size() ; c++) {
		if (atomList.atomAt(c).name().equals("CA") && (atomList.atomAt(c).residueNumber()>=begins) && 
		(atomList.atomAt(c).residueNumber()<=ends)) {
			atoms[cc] = atomList.atomAt(c);
			cc++;
		}
	}
	comment = "localRG";
    }
 
    public void evaluateAtoms() {}
    
    public double evaluate() {
	if (! on) return 0.0;
 	double energy = 0;
	for (int c=0;c<distances.length ; c++)
		distances[c] = 1e10;
 	Atom atom;
 	double dis;
 	for (int c=0 ; c<atomList.size() ; c++) {
 		atom = atomList.atomAt(c);
 		if (!atom.isHydrogen && ((atom.residueNumber()>ends) || (atom.residueNumber()<begins)))
 		for (int d=0 ; d<atoms.length ; d++) {
 			dis = (atoms[d].x()-atom.x())*(atoms[d].x()-atom.x()) +
 				(atoms[d].y()-atom.y())*(atoms[d].y()-atom.y()) +
 				(atoms[d].z()-atom.z())*(atoms[d].z()-atom.z());
 			if (dis < distances[d])
 				distances[d] = dis;
 		}
 	}
 	for (int d=0 ; d<atoms.length ; d++) 	 
 		energy += Math.sqrt(distances[d]);
 			
    return energy;
 	} 
 	
}