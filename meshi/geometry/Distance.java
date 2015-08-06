package meshi.geometry;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.AtomPair;
import meshi.util.Attributable;
import meshi.util.AttributesRack;
import meshi.util.MeshiAttribute;
import meshi.util.filters.Filter;
import meshi.util.formats.Fdouble;
/** 
 * The distance between two {@link meshi.molecularElements.Atom Atoms}<br>.
 *
 * Almost any measurable feature of a molecule is related to distances
 * between pairs of {@link meshi.molecularElements.Atom Atoms}.
 * Thus, The calculation of distances (and their inverse and derivatives)
 * are typically a computational bottleneck in computational structural 
 * biology applications. In all applications that we are aware of (Please,
 * enlighten us if you know better) distance calculation is done as part 
 * of the procedures that use it (say, as part of the van-der-Waals energy 
 * calculation). As a result the distance between two atoms may be calculated 
 * more then once. For example the distance between two atoms may be 
 * calculated both during angle and torsion angle energies calculations.
 * In Meshi We tried to consentrate all distance related issues in a few 
 * classes: this one, its subclasses and the closely connected 
 * class  {@link meshi.geometry.DistanceMatrix DistanceMatrix}.
 * <p>
 * <b> Possible pitfalls </b><br>
 * <ol>
 *  <li>The object variables dx, dy, dz, distance, distance2, invDistance, 
 *      dDistanceDx, dDistanceDy and dDistanceDz should have been private and 
 *      accessed through "get methods". Energy functions though, use the 
 *      values of these variables intensively and a considerable computational
 *       gain is achieved granting them public accessability (protected for 
 *       dx, dy and dz) and removing the function call overhead. Thus, 
 *       changing the values of these variables by other classes is 
 *       possible (from the compiler point of view) but not very much 
 *       recommended.
 *  <li> The distance object <b>"is not aware" </b> of changes in the 
 *       atoms coordinates. The values stored in this object are correct 
 *       only after the {@link meshi.geometry.Distance#update update}
 *       method is <u>explicitly</u> called.  
 **/
public class Distance implements Comparable, Attributable{
    //---------------------------- fields ------------------------------
    public static final double INFINITE_DISTANCE = Double.MAX_VALUE;
    /**
     * This object represent the distance between atom1 and atom2.
     **/
    protected Atom atom1, atom2;
    protected int atom1Number;
    protected int atom2Number;
    protected final int atom1Number() { return atom1Number;}
    protected final int atom2Number() { return atom2Number;}
    public static final Filter filter = new IsDistance();
	
    /**
     * The inverse of the distance between atom1 & atom2.         <br>
     * Not that the public accessability of of this variable
     * improves computational efficiency but opens a wide door 
     * for bugs.
     **/
    private double invDistance;
	
    /**
     * The derivative of the distance between atom1 & atom2 by 
     * the X coordinate of atom1. The derivative by the X coordinate of 
     * atom2 is this value multiplied by -1.                      <br>
     * Not that the public accessability of of this variable
     * improves computational efficiency but opens a wide door 
     * for bugs.
     **/
    private double dDistanceDx;
	
    /**
     * The derivative of the distance between atom1 & atom2 by 
     * the Y coordinate of atom1. The derivative by the Y coordinate of 
     * atom2 is this value multiplied by -1.                      <br>
     * Not that the public accessability of of this variable
     * improves computational efficiency but opens a wide door 
     * for bugs.
     **/
    private double dDistanceDy;
	
    /**
     * The derivative of the distance between atom1 & atom2 by 
     * the Z coordinate of atom1. The derivative by the Z coordinate of 
     * atom2 is this value multiplied by -1.                      <br>
     * Not that the public accessability of of this variable
     * improves computational efficiency but opens a wide door 
     * for bugs.
     **/
    private double dDistanceDz;
	
    /**
     * The square of the distance between atom1 & atom2  <br>
     * Not that the public accessability of of this variable
     * improves computational efficiency but opens a wide door 
     * for bugs.
     **/
    protected double distance2;
	
    /**
     * The distance between atom1 & atom2  <br>
     * Not that the public accessability of of this variable
     * improves computational efficiency but opens a wide door 
     * for bugs.
     **/
    protected double distance;
    /**
     * The distance between atom1 & atom2 on the previous step <br>
     **/
    private double distancePrev = INFINITE_DISTANCE;
	
    /**
     * atom1.x - atom2.x
     * Not that the protected accessability of of this variable
     * improves computational efficiency but opens a wide door 
     * for bugs.
     **/
    private double dx;
    /**
     * atom1.y - atom2.y
     * Not that the protected accessability of of this variable
     * improves computational efficiency but opens a wide door 
     * for bugs.
     **/
    private double dy;
    /**
     * atom1.z - atom2.z
     * Not that the protected accessability of of this variable
     * improves computational efficiency but opens a wide door 
     * for bugs.
     **/
    private double dz;
	
    /**
     * The coordinates of atom1.
     **/
    protected Coordinates coor1;
	
    /**
     * The coordinates of atom2.
     **/
    protected Coordinates coor2;
    //protected double rMax = DistanceMatrix.rMax();    
    //protected double rMax2 = DistanceMatrix.rMax2();    
    //protected double rMaxPlusBuffer2 = DistanceMatrix.rMaxPlusBuffer2();      
    private int debugFlag = -1;
    protected boolean bonded = false;
    public boolean bonded() {return bonded;}
    protected void setBonded() {
	bonded = true;
    }
    protected boolean mirror;
    protected boolean mirror() {return mirror;}
    public final boolean frozen;
    protected int largeType, smallType;
    public int largeType()  {return largeType;}
    public int smallType()  {return smallType;}
    private AtomList atoms;
    private AttributesRack attributes = new AttributesRack();
    public final void addAttribute(MeshiAttribute attribute){ attributes.addAttribute(attribute); }
    public final MeshiAttribute getAttribute(int key){ return attributes.getAttribute(key); }
    
    public final AtomList atoms() {return atoms;}
    //-------------------------- Constructors --------------------------
    public Distance() {
	distance = INFINITE_DISTANCE;
	bonded = false;
	debugFlag =1;
	frozen = true;
    }
    /**
     * Instantiate and update a Distance Object.
     **/
     public Distance(AtomPair atomPair) {
	this(atomPair.atom1(),atomPair.atom2());
    }
	
    /**
     * Instantiate and update a Distance Object.
     **/
    public  Distance(Atom atom1, Atom atom2) {
	this.atom1 = atom1;
	this.atom2 = atom2;
	atom1Number = atom1.number();
	atom2Number = atom2.number();
	coor1 = atom1.coordinates();
	coor2 = atom2.coordinates();    
	update();        
	bonded = false;
	if (atom1.frozen() & atom2.frozen()) frozen = true;
	else frozen = false;
	mirror = false;
       if (atom1.iType() > atom2.iType()) {
                largeType = atom1.iType();
                smallType = atom2.iType();
        }
        else {
                largeType = atom2.iType();
                smallType = atom1.iType();
        }
        atoms = new AtomList(AtomPair.ATOM_PAIR_CAPACITY);
        atoms.fastAdd(atom1);
        atoms.fastAdd(atom2);
    }
	
    /**
     * Instantiate and update a Distance Object.
     * This constructor should be used with caution. It helps to save time in DistanceMatrix update. 
     * We see no reason to use it anywhere else.  
     **/
    protected  Distance(Atom atom1, Atom atom2, double realD2, double dx, double dy, double dz,double rMax2) {
	this.atom1 = atom1;
	this.atom2 = atom2;
	atom1Number = atom1.number();
	atom2Number = atom2.number();
	coor1 = atom1.coordinates();
	coor2 = atom2.coordinates();	         
	distance2 = realD2;    
	if (distance2 < rMax2) {
	    this.dx = dx;
	    this.dy = dy;
	    this.dz = dz;
	    distance = Math.sqrt(distance2);
	    invDistance = 1/distance;
	    dDistanceDx = dx*invDistance;
	    dDistanceDy = dy*invDistance;
	    dDistanceDz = dz*invDistance;
	}
	else {
	    distance = INFINITE_DISTANCE;
	    debugFlag = 2;
	}
	bonded = false;
	if (atom1.frozen() & atom2.frozen()) frozen = true;
	else frozen = false;
	mirror = false;
       if (atom1.iType() > atom2.iType()) {
                largeType = atom1.iType();
                smallType = atom2.iType();
        }
        else {
                largeType = atom2.iType();
                smallType = atom1.iType();
        }
        atoms = new AtomList(AtomPair.ATOM_PAIR_CAPACITY);
        atoms.fastAdd(atom1);
        atoms.fastAdd(atom2);
    }                        
	
	
	
    //------------------------------- methods ------------------------------
	
    public boolean equals(Object obj) {
	    Distance other = (Distance) obj;
	   boolean y = ((atom1.equals(other.atom1) && atom2.equals(other.atom2)) ||
                    (atom1.equals(other.atom2) && atom2.equals(other.atom1)));
        return y;
   }

      public boolean amiEquals(Object obj) {
	    Distance other = (Distance) obj;
	   boolean y = ((atom1.number() == other.atom1Number && atom2.number() == other.atom2Number) ||
                    (atom1.number() == other.atom2Number && atom2.number() == other.atom1Number));
        return y;
    }
    
    //method equals for NonBondedList testing 
    //we need the spesial checking (instead ==) because attributes of distances are not equal any more
    public boolean testEquals(Distance other) {
        return (atom1.equals(other.atom1) && atom2.equals(other.atom2) &&
                        distance == other.distance);
    }

    /**
     * The distance values.
     **/ 
    public boolean update() {
	return update(1000000,1000000);
    }
    /**
     * The distance values.
     * This method should be used with caution. We see no point in using it outside DistanceMatrix. 
     **/
 
    protected boolean update(double rMax2, double rMaxPlusBuffer2) { 
	distancePrev = this.distance;	
        try {
		dx  = coor1.x[0] - coor2.x[0];	
		dy  = coor1.y[0] - coor2.y[0];	
		dz  = coor1.z[0] - coor2.z[0];
        }
	catch (RuntimeException ex) { System.out.println("xxxxx "+this); throw ex;}
	distance2 = dx*dx + dy*dy + dz*dz;              
	if ((distance2 < rMax2) || bonded) { 
	    distance = Math.sqrt(distance2);	    
	    invDistance = 1/distance;
	    dDistanceDx = dx*invDistance;
	    dDistanceDy = dy*invDistance;
	    dDistanceDz = dz*invDistance;
	    return true;
	}
	distance = INFINITE_DISTANCE;	
	if (distance2 < rMaxPlusBuffer2) {                    
	    debugFlag = 4;
	    return true;
	}
	return false;
    }
	
    /**
     * Get the distance between atom1 and atom2.
     **/
    public  double distance() { 
	return distance;
    }
	
    /**
     * Get the square of the distance between atom1 and atom2.
     **/
    public  double distance2() { return distance2;}
	
    /**
     * Get the inverse of the distance between atom1 and atom2.
     **/
    public  double invDistance() {return invDistance;}
    
    /**
     * Get the previous of the distance between atom1 and atom2.
     **/
    public  double distancePrev() {return distancePrev;}
	
    /**
     * The derivative of the distance between atom1 & atom2 by 
     * the X coordinate of atom1. The derivative by the X coordinate of 
     * atom2 is this value multiplied by -1.                     
     **/    
    public  double dDistanceDx() {return dDistanceDx;}
	
    /**
     * The derivative of the distance between atom1 & atom2 by 
     * the Y coordinate of atom1. The derivative by the Y coordinate of 
     * atom2 is this value multiplied by -1.                     
     **/    
    public  double dDistanceDy() {return dDistanceDy;}
	
    /**
     * The derivative of the distance between atom1 & atom2 by 
     * the Z coordinate of atom1. The derivative by the Z coordinate of 
     * atom2 is this value multiplied by -1.                     
     **/    
    public  double dDistanceDz() {return dDistanceDz;}
	
    /**
     * Get atom1.
     **/
    public  Atom atom1() { return atom1; }
	
    /**
     * Get atom2.
     **/
    public  Atom atom2() { return atom2; }
	
    public  double dx() {return dx;}
	
    public  double dy() {return dy;}
	
    public  double dz() {return dz;}    
	
    public String toString() {
	Fdouble dformat = Fdouble.STANDARD;
	return "Distance "+atom1.number()+" "+atom2.number()+" "+
	    dformat.f(distance)+" "+
	    dformat.f(dDistanceDx)+" "+ 
	    dformat.f(dDistanceDy)+" "+ 
	    dformat.f(dDistanceDz);
    }
    public String comment() { return "Distance";}

    public boolean isDistanceOf(Atom toCueck){
        return atom1.equals(toCueck ) | atom2 .equals(toCueck );
    }

	
    public int compareTo(Object obj) {
	Distance other = (Distance) obj;
	if  (atom2Number > other.atom2Number) return 1;
	if (atom2Number < other.atom2Number) return -1;
	return 0;
    }
    private static class IsDistance implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof Distance);
	}
    } 
}

