package meshi.geometry.hydrogenBond;

import java.text.DecimalFormat;

import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.Activable;


/** 
 * A super class with the most general treatment of hydrogen bond. 
 * A specific implementation should put content to these abstract methods.
 * The constructor gets a distance matrix object, and an AtomList object of 
 * all the atoms participating in the bond. The only condition on the list is
 * that the donor and acceptor must be its first two elements (their relative
 * order is not important).
 *
 * IPORTANT IMPORTANT IMPORTANT
 * NOTE: The distance matrix must be in an updated state before any of this class's
 * methods are called. If it is not up to date, the results would be meaningless!!!
 * This class does not update the distance matrix on its own at any time!!
 *
 **/
public abstract class AbstractHydrogenBond implements Activable {

    protected AtomList atomList=null;
    protected DistanceMatrix dm=null; 
    protected double hbVal = 0.0;
    protected boolean active;

    /** 
    * derivatives[i][0] is the derivative of the HB value relative to the X coordinate of atom i.
    * derivatives[i][1] is the derivative of the HB value relative to the Y coordinate of atom i.
    * derivatives[i][2] is the derivative of the HB value relative to the Z coordinate of atom i.
    **/
    protected double[][] derivatives = null;

    public AbstractHydrogenBond() {}

    public AbstractHydrogenBond(DistanceMatrix dm, AtomList atomList) {
    	this.dm = dm;
    	this.atomList = atomList;
    	derivatives = new double[atomList.size()][3];
    	addObjectToActivableListInAtoms();
    	updateActivity();
    }

	public void addObjectToActivableListInAtoms() {
		for (int c=0; c<atomList.size() ; c++)
			atomList.atomAt(c).addActivable(this);
	}
    
	public void updateActivity() {
		active = false;
		for (int c=0; c<atomList.size() ; c++)
			if (atomList.atomAt(c).active())
				active = true;
	}
	
	public boolean active() {
		return active;
	}
    
    /** 
      * This update method should be used if a change to the atoms coordinates has occured.
      * Only after this update, the getter methods can produce correct values.
      * This class assume the distance matrix is updated. 
      **/
    protected abstract void updateHBvalueAndDerivatives();
    

    /** 
      * FOLLOWING the activation of "update" this method can be activated, 
      * where it will apply the forces to the atoms participating in the HB.
      * The forces are: factor*(minus derivative). This class is useful in energy terms that
      * make use of the hydrogen bond. 
      **/
	public void applyForcesToAtoms(double factor) {
		Atom atom;
		
		for (int c=0 ; c<atomList.size() ; c++) {
			atom = atomList.atomAt(c);
			if (! atom.frozen())
			atomList.atomAt(c).addToFx(-factor*derivatives[c][0]);
			atomList.atomAt(c).addToFy(-factor*derivatives[c][1]);
			atomList.atomAt(c).addToFz(-factor*derivatives[c][2]);
		}
	}
	
	public double hbVal() {return hbVal;}

	public Atom getFirstPolar() {return atomList.atomAt(0);}	
	public Atom getSecondPolar() {return atomList.atomAt(1);}	

	public String toString() {
		DecimalFormat format = new DecimalFormat("0.##");
		return "HB: " + atomList.atomAt(0).residue() + " " + atomList.atomAt(0).name() + " --- " + 
		atomList.atomAt(1).residue() + " " + atomList.atomAt(1).name() + "   hbVal: " + format.format(hbVal()); 
	}
}


  
