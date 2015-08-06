package meshi.geometry;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomPair;
import meshi.util.Updateable;
import meshi.util.UpdateableException;
/**
 *
 *--------------------------- Angle -----------------------------------
 **/
	
public class Angle  implements Updateable {
    public final Atom atom1;
    public final Atom atom2;
    public final Atom atom3;
    protected String angleName = "";
    protected int angleCode = -1;
    protected int angleResNum = -9999;
    protected String angleResName = "";    
    protected AtomPair atomPair1, atomPair2;
    public final Distance DISTANCE1, DISTANCE2;
    private double angle,fakeAngle;
    private double dangleDx1, dangleDy1, dangleDz1;
    private double dangleDx2, dangleDy2, dangleDz2;
    private double dangleDx3, dangleDy3, dangleDz3;
    private double cosAngle;
    private double sinAngle;
    private double dFakeAngledAngle = 1.0;
    protected DistanceMatrix distanceMatrix = null;
    private boolean useFakeAngle = false;
    private int numberOfUpdates = 0;

    /* Once upon a time there were BUG GENERATING constructors that did not use distance matrix.
       After they made enough damage they were removed */


    public Angle(Atom atom1, Atom atom2, Atom atom3, DistanceMatrix distanceMatrix, boolean useFakeAngle) {
	this(new AtomPair(atom1, atom2), new AtomPair(atom2, atom3), distanceMatrix, useFakeAngle); 
    } 
    
    public Angle(AtomPair atomPair1, AtomPair atomPair2, 
		 DistanceMatrix distanceMatrix) {
	this(atomPair1, atomPair2, distanceMatrix, false);
    }

    /**
     * This constructor should be used only by DisposableAngle.
     **/
    protected Angle(Atom atom1, Atom atom2, Atom atom3, Distance distance1, Distance distance2) {
	this.useFakeAngle = false;
	this.distanceMatrix = null;
	this.atom1 = atom1;
	this.atom2 = atom2;
	this.atom3 = atom3;

	DISTANCE1 = distance1;
	DISTANCE2 = distance2;
	assignName();
        update();
    }

    public Angle(AtomPair atomPair1, AtomPair atomPair2, 
		 DistanceMatrix distanceMatrix,boolean useFakeAngle) {
        this.useFakeAngle = useFakeAngle;
	this.atomPair1 = atomPair1;
	this.atomPair2 = atomPair2;
	this.distanceMatrix = distanceMatrix;
	Distance distance1, distance2;
	Distance distance1mirror, distance2mirror;

	distance1 = distanceMatrix.distance(atomPair1.atom1(), 
					    atomPair1.atom2());
	distance2 = distanceMatrix.distance(atomPair2.atom1(), 
					    atomPair2.atom2());
	   
	distance1mirror = distanceMatrix.distance(atomPair1.atom2(), 
						  atomPair1.atom1());
	distance2mirror = distanceMatrix.distance(atomPair2.atom2(), 
						  atomPair2.atom1());

    if (((atomPair1.atom1() == atomPair2.atom1()) && (atomPair1.atom2() == atomPair2.atom2())) ||
	    ((atomPair1.atom1() == atomPair2.atom2()) && (atomPair1.atom2() == atomPair2.atom1())))
	   throw new RuntimeException("Cannot create an angle from "+
					atomPair1+
					" and "+atomPair2);

	if (atomPair1.atom1() == atomPair2.atom1()) { 
	    //  2-1
	    //    1-2
	    atom1 = atomPair1.atom2();
	    atom2 = atomPair1.atom1();
	    atom3 = atomPair2.atom2();
	    DISTANCE1 = distance1mirror;
	    DISTANCE2 = distance2mirror;
	}
	else if (atomPair1.atom1() == atomPair2.atom2()) { 
	    // 2-1
            //   2-1 
	    atom1 = atomPair1.atom2();
	    atom2 = atomPair1.atom1();
	    atom3 = atomPair2.atom1();
	    DISTANCE1 = distance1mirror;
	    DISTANCE2 = distance2;
	}
	else if (atomPair1.atom2() == atomPair2.atom2()) { 
	    // 1-2
	    //   2-1
	    atom1 = atomPair1.atom1();
	    atom2 = atomPair1.atom2();
	    atom3 = atomPair2.atom1();
	    DISTANCE1 = distance1;
	    DISTANCE2 = distance2;
	}
	else if (atomPair1.atom2() == atomPair2.atom1()) { 
	    // 1-2
	    //   1-2
	    atom1 = atomPair1.atom1();
	    atom2 = atomPair1.atom2();
	    atom3 = atomPair2.atom2();
	    DISTANCE1 = distance1;
	    DISTANCE2 = distance2mirror;
	}
	else throw new RuntimeException("Cannot create an angle from "+
					atomPair1+
					" and "+atomPair2);
	assignName();
    update();
    }

    public Atom atom1() {return atom1;}
    public Atom atom2() {return atom2;}
    public Atom atom3() {return atom3;}
    public AtomPair atomPair1() {return atomPair1;}
    public AtomPair atomPair2() {return atomPair2;}

    public AtomPair sharedAtomPair(Angle other) {
	if (atomPair1.equals(other.atomPair1)) return atomPair1;
	if (atomPair1.equals(other.atomPair2)) return atomPair1;
	if (atomPair2.equals(other.atomPair1)) return atomPair2;
	if (atomPair2.equals(other.atomPair2))  return atomPair2;
	return null;
    }
    public boolean proper(Angle other) {
	AtomPair shared = sharedAtomPair(other);
	if (shared == null) return false;
	if (atom2 == other.atom2) return false;
	return true;
    }
    public AtomPair hinge(Angle other) {
	AtomPair out = sharedAtomPair(other);
	if (out == null) return null;
	if ((atom2 == other.atom1()) |
	    (atom2 == other.atom3())) return out;
	return null;
    }
    public String toString (){ 
	return "Angle: "+ angleName + 
	    "\n\tatom1   :"+atom1()+" "+atom1().getClass()+
	    "\n\tatom2   :"+atom2()+" "+atom2().getClass()+
	    "\n\tatom3   :"+atom3()+" "+atom3().getClass()+
	    "\n\tvalue   :"+rad2deg(angle());
    }
    public double angle() {
        if(useFakeAngle)
            return fakeAngle;
        return angle;
    }
    public double dangleDx1() { return dangleDx1;}
    public double dangleDy1() { return dangleDy1;}
    public double dangleDz1() { return dangleDz1;}
    public double dangleDx2() { return dangleDx2;}
    public double dangleDy2() { return dangleDy2;}
    public double dangleDz2() { return dangleDz2;}
    public double dangleDx3() { return dangleDx3;}
    public double dangleDy3() { return dangleDy3;}
    public double dangleDz3() { return dangleDz3;}
    public double cosAngle() {  return cosAngle;}
    public double sinAngle() { return sinAngle;}
    public double dFakeAngledAngle() { return dFakeAngledAngle;}

    public void update(int numberOfUpdates) throws UpdateableException{
	if (numberOfUpdates == this.numberOfUpdates+1) {
	    update();
	    this.numberOfUpdates++;
	}
	else if (numberOfUpdates != this.numberOfUpdates) 
	    throw new RuntimeException("Something weird with Angle.update(int numberOfUpdates)\n"+
				       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
    }

   protected void update() {
	double dd1Dx = DISTANCE1.dDistanceDx();
	double dd1Dy = DISTANCE1.dDistanceDy();
	double dd1Dz = DISTANCE1.dDistanceDz();
	double invD1 = DISTANCE1.invDistance();
	double dd2Dx = -1*DISTANCE2.dDistanceDx();
	double dd2Dy = -1*DISTANCE2.dDistanceDy();
	double dd2Dz = -1*DISTANCE2.dDistanceDz();
	double invD2 = DISTANCE2.invDistance();
	
	cosAngle = (dd1Dx*dd2Dx + dd1Dy*dd2Dy + dd1Dz*dd2Dz);
        // It is assumed that these values will never occur in normal run. 
	// They may occur in some weird initialization situation.
	if (cosAngle >= 1.0) {
		angle = Math.PI;
		dangleDx1 = 0;
		dangleDy1 = 0;
		dangleDz1 = 0;
		dangleDx2 = 0;
		dangleDy2 = 0;
		dangleDz2 = 0;
		dangleDx3 = 0;
		dangleDy3 = 0;
		dangleDz3 = 0;		
	}		
	else if (cosAngle <= -1.0) {
		angle = 0;
		dangleDx1 = 0;
		dangleDy1 = 0;
		dangleDz1 = 0;
		dangleDx2 = 0;
		dangleDy2 = 0;
		dangleDz2 = 0;
		dangleDx3 = 0;
		dangleDy3 = 0;
		dangleDz3 = 0;	
	}	
	else {
            angle = ArcCos.acos(-1*cosAngle);//runs at least four times faster
	    //            angle = Math.acos(-1*cosAngle);
            //                 {angle                                          if angle <= 170 
            //fakeAngle(angle)={-0.008*angle^3 +4.04*angle^2 -679*angle+38148  if 170<angle<= 175
            //                 {173                                            if angle > 175 
            fakeAngle = angle;
            dFakeAngledAngle = 1.0;
                   
            if (useFakeAngle
                && (angle > deg2rad(170)) && (angle <= deg2rad(175))){
                double angleDeg = rad2deg(angle);
		fakeAngle = deg2rad(-0.008*angleDeg*angleDeg*angleDeg +4.04*angleDeg*angleDeg -679*angleDeg+38148);
                cosAngle = Math.cos(fakeAngle);
		dFakeAngledAngle = deg2rad(3*(-0.008)*angleDeg*angleDeg+2*4.04*angleDeg-679);//why do we need this - ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
	    }
	    else {
		if (useFakeAngle && angle > deg2rad(175)) {
                    fakeAngle = deg2rad(173);
                    cosAngle = Math.cos(fakeAngle);
		    dFakeAngledAngle = 0.0;
		}
	    }
            sinAngle = Math.sin(fakeAngle);
            double factor = invD1/sinAngle;
            dangleDx1 = factor * (dd2Dx - cosAngle * dd1Dx);
            dangleDy1 = factor * (dd2Dy - cosAngle * dd1Dy);
            dangleDz1 = factor * (dd2Dz - cosAngle * dd1Dz);
            factor = invD2/sinAngle;
            dangleDx3 = factor * (cosAngle * dd2Dx - dd1Dx);
            dangleDy3 = factor * (cosAngle * dd2Dy - dd1Dy);
            dangleDz3 = factor * (cosAngle * dd2Dz - dd1Dz);
            dangleDx2 = -1*(dangleDx1 + dangleDx3);
            dangleDy2 = -1*(dangleDy1 + dangleDy3);
            dangleDz2 = -1*(dangleDz1 + dangleDz3);
            
	}
    }

    public boolean frozen() {
	return (atom1.frozen() & atom2.frozen() & atom3.frozen());
    }
    public static double deg2rad(double ang) {
	return ang*Math.PI/180;
    }
    public static double rad2deg(double ang) {
	return ang*180/Math.PI;
    }
    
    
    public String getAngleName() {return angleName;}
    public int getAngleCode() {return angleCode;}
    public int getAngleResNum() {return angleResNum;}
    public String getAngleResName() {return angleResName;}
    /**
     * Assigns a meaningful name to an angle.
     * This method assigns (if possible):
     *  - name to the angle object.
     *  - number corresponding to this name.
     *  - the number in the chain of the center atom residue to which this amgle belong
     *  - the name of the center atom residue to which this angle belong 
     *
     *The names to numbers conversion is:
     *Three consecutive Ca's - 0
     **/    
    protected void assignName() {
    	int resNum = atom2.residueNumber();
    	if ((atom1.name().compareTo("CA") == 0) &&                 // CA3
    		(atom2.name().compareTo("CA") == 0) && 
    		(atom3.name().compareTo("CA") == 0) && 
    		(((atom3.residueNumber() == (resNum+1)) && (atom1.residueNumber() == (resNum-1))) ||
    		 ((atom1.residueNumber() == (resNum+1)) && (atom3.residueNumber() == (resNum-1))))){
     			angleName = "CA3";
     			angleCode = 0;
     			angleResNum = resNum;
     			angleResName = atom2.residueName();
     	}
     }    
}
