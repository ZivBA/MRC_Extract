package meshi.util;
import java.util.Iterator;

import meshi.geometry.Coordinates;
import meshi.molecularElements.Atom;
import meshi.molecularElements.Protein;
import meshi.sequences.AtomAlignment;
import meshi.sequences.ResidueAlignment;
import meshi.util.overlap.Overlap;


/**
 * A (hopefully) comfortable handle to the complexity of the Kabsch structural overlap algorithm[1,2]. 
 * The algorithm itself is implemented in meshi.overlap.Overlap.
 * 1. "A solution for the best rotation to relate two sets of vectors". By Wolfgang Kabsch, (1976) Acta Cryst. A 32: 922-923<BR>
 * 2. "A discussion of the solution for the best rotation to relate two sets<BR>
 *    of vectors". By Wolfgang Kabsch, Acta Cryst. (1978). A34, 827-828.
 * This class was rather drastically modified by Chen in version 1.45.
 **/
public class Rms implements KeyWords{

    private static boolean debug = true;
    private boolean alive = true;
    //data members
    private double[][] rotateMatrix;
    private double rms;
    private Coordinates centerOfMass0 = new Coordinates();
    private Coordinates centerOfMass1 = new Coordinates();
        

    //constructors

    public Rms(ResidueAlignment residueAlignment) {
	this(new AtomAlignment( residueAlignment));
    }
	
    public Rms(AtomAlignment atomAlignment) {
	if (atomAlignment.hasGaps()) throw new RuntimeException("Cannot calculate RMS for AtomAlignment with gaps");
	int size = atomAlignment.size();
	double[][] coor0 = new double[3][size];
	double[][] coor1 = new double[3][size];
	String comment0 = atomAlignment.columnAt(0).cell(0).comment;
	String comment1 = atomAlignment.columnAt(0).cell(1).comment;
	double centerOfMassX1 = 0, centerOfMassY1 = 0, centerOfMassZ1 = 0;
	double centerOfMassX0 = 0, centerOfMassY0 = 0, centerOfMassZ0 = 0;
	try {
	    for (int i = 0; i < size; i++){
		Atom atom0 = atomAlignment.atomAt(i,0);
		Atom atom1 = atomAlignment.atomAt(i,1);
		coor0[0][i] = atom0.x();
		coor0[1][i] = atom0.y();
		coor0[2][i] = atom0.z();
		coor1[0][i] = atom1.x();
		coor1[1][i] = atom1.y();
		coor1[2][i] = atom1.z();
		
		centerOfMassX0 += atom0.x();
		centerOfMassY0 += atom0.y();
		centerOfMassZ0 += atom0.z();
		centerOfMassX1 += atom1.x();
		centerOfMassY1 += atom1.y();
		centerOfMassZ1 += atom1.z();
	    }
	    centerOfMassX0 /= size;
	    centerOfMassY0 /= size;
	    centerOfMassZ0 /= size;
	    centerOfMassX1 /= size;
	    centerOfMassY1 /= size;
	    centerOfMassZ1 /= size;

	    centerOfMass0.setX(centerOfMassX0);
	    centerOfMass0.setY(centerOfMassY0);
	    centerOfMass0.setZ(centerOfMassZ0);
	    centerOfMass1.setX(centerOfMassX1);
	    centerOfMass1.setY(centerOfMassY1);
	    centerOfMass1.setZ(centerOfMassZ1);
	    
	    Overlap overlap = new Overlap(coor0, coor1, size, comment0, comment1);
	    rms = overlap.rms();
	    rotateMatrix = overlap.rotationMatrix();
	}
	catch (Exception e) {
	    throw new MeshiException("Rms Error:\n"+
				  "comparing "+comment0+"\n"+
				  "with\n"+comment1+"\n"+
				  "reproted problem"+
				  e);
	}
    }
	

    private Coordinates centerOfMass0() { return centerOfMass0;}
    private Coordinates centerOfMass1() { return centerOfMass1;}

    public static double superimpose(Protein protein0, Protein protein1, ResidueAlignment residueAlignment){
	AtomAlignment atomAlignment = new AtomAlignment(residueAlignment);
	 Rms rms = new Rms(atomAlignment);
	 Coordinates centerOfMass0 = rms.centerOfMass0();
	 Coordinates centerOfMass1 = rms.centerOfMass1();
	 for (Iterator atoms = protein1.atoms().iterator(); atoms.hasNext();) {
	     Atom atom = (Atom) atoms.next();
	     atom.addToX(0-centerOfMass1.x());
	     atom.addToY(0-centerOfMass1.y());
	     atom.addToZ(0-centerOfMass1.z());
	 }
	 double[][] rotateMatrix = rms.getMatrix();
	 for (Iterator atoms = protein1.atoms().iterator(); atoms.hasNext();) {
	     Atom atom = (Atom) atoms.next();
	     double x = atom.x();
	     double y = atom.y();
	     double z = atom.z();
	     atom.setX(rotateMatrix[0][0]*x + rotateMatrix[0][1]*y + rotateMatrix[0][2]*z);
	     atom.setY(rotateMatrix[1][0]*x + rotateMatrix[1][1]*y + rotateMatrix[1][2]*z);
	     atom.setZ(rotateMatrix[2][0]*x + rotateMatrix[2][1]*y + rotateMatrix[2][2]*z);
	 }
	 for (Iterator atoms = protein1.atoms().iterator(); atoms.hasNext();) {
	     Atom atom = (Atom) atoms.next();
	     atom.addToX(centerOfMass0.x());
	     atom.addToY(centerOfMass0.y());
	     atom.addToZ(centerOfMass0.z());
	 }
	 return rms.getRms();
    }
	 	
	




    //methods

    public double[][] getMatrix(){
	if (alive) return rotateMatrix;
	else return null;
    }

  
    public double getRms(){
    	if (alive) return rms;
	else return -1;
    }

    public String toString() {
	return (new Double(getRms())).toString();
    }

    public Coordinates getCenterOfMass0() {return centerOfMass0;}
    public Coordinates getCenterOfMass1() {return centerOfMass1;}

}
