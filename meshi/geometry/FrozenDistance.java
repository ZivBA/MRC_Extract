package meshi.geometry;
import meshi.molecularElements.Atom;
public class FrozenDistance extends Distance{

//     public FrozenDistance(AtomPair atomPair) {
// 	this(atomPair.atom1(),atomPair.atom2());
//      	super.update();
//     	super.update(1000000,10000000);
//     }
//     protected  FrozenDistance(Atom atom1, Atom atom2, double realD2, double realDx, double realDy, double realDz,double rMax2) {
// 	super (atom1, atom2, realD2, realDx, realDy, realDz, rMax2);
// 	if ((! atom1.frozen()) | (!atom2.frozen())) 
// 		throw new RuntimeException("Only frozen atoms use frozen distance");
//     	super.update(1000000,10000000);
//      }	

//    public FrozenDistance(Atom atom1, Atom atom2) {
//     	super(atom1,atom2);
// 	if ((! atom1.frozen()) | (!atom2.frozen())) 
// 		throw new RuntimeException("Only frozen atoms use frozen distance");
//     	super.update(1000000,10000000);
//     }
   public FrozenDistance(Atom atom1, Atom atom2,double rMax2, double rMaxPlusBuffer2) {
    	super(atom1,atom2);
	if ((! atom1.frozen()) | (!atom2.frozen())) 
		throw new RuntimeException("Only frozen atoms use frozen distance");
    	super.update(rMax2,rMaxPlusBuffer2);
    }
	
    //------------------------------- methods ------------------------------

    /**
     * The distance values.
     **/ 
    public boolean update() {
    	return true;
    }

    protected boolean update(double rMax2, double rMaxPlusBuffer2) { 
	return true;
    }
	
		
	
	
    public String comment() { return "FrozenDistance";}
}

