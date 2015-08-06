package meshi.geometry;
import meshi.molecularElements.Atom;
public class DistanceMirror  extends Distance {
    public final Distance source;
    public  DistanceMirror(Distance source) {
	super();
	this.source = source;
      atom2Number = source.atom1.number();
      atom1Number = source.atom2.number();
      if (atom2Number < atom1Number) throw new RuntimeException("weird distance mirror \n"+
		      	"atom2Number = "+atom2Number+"  tom1Number = "+atom1Number+"\n"+
				atom1()+"\n"+atom2());
      mirror = true;
    }    
    public final Atom atom1() {return source.atom1();}
    public final Atom atom2() {return source.atom2();}
    public final double distance() {return source.distance();}
    public final double distance2() { return source.distance2();}
    public final double invDistance() { return source.invDistance();}
    public final double dDistanceDx() { return 0 - source.dDistanceDx();}
    public final double dDistanceDy() { return 0 - source.dDistanceDy();}
    public final double dDistanceDz() { return 0 - source.dDistanceDz();}
    public final double dx() { return 0 - source.dx();}
    public final double dy() { return 0 - source.dy();}
    public final double dz() { return 0 - source.dz();}

   public String toString() {
	return "Distance mirror of "+source.comment()+" "+
	    source.atom1().number()+" "+
	    source.atom2().number()+" "+
	    distance()+" "+
	    dDistanceDx()+" "+
	    dDistanceDy()+" "+
	    dDistanceDz()+" ";
    }
    public String comment() {return "Distance mirror";}
    public boolean update() {
	    throw new RuntimeException("cannot update distanceMirror "+this);
    }
}

