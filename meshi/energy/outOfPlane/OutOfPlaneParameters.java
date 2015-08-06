package meshi.energy.outOfPlane;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.geometry.Angle;
import meshi.molecularElements.Atom;
import meshi.util.filters.Filter;
public class OutOfPlaneParameters extends Parameters implements Comparable {
    public final double target, force, force2;
    public final int type1, type2, type3, type4;
    public OutOfPlaneParameters() {
	this(-1,-1,-1,-1,-9999.9,-9999.9); 
    }
    public OutOfPlaneParameters(int type1, int type2, int type3, int type4) {
	this(type1,type2,type3,type4,-9999.9,-9999.9); 
    }
    public OutOfPlaneParameters(String line) {
	this(new StringTokenizer(line));
    }
    private OutOfPlaneParameters(StringTokenizer line) {
	this(Atom.type(line.nextToken()), 
	     Atom.type(line.nextToken()), 
	     Atom.type(line.nextToken()), 
	     Atom.type(line.nextToken()), 
	     toDouble(line.nextToken()), toDouble(line.nextToken()));
    }

    public OutOfPlaneParameters(int type1, int type2,
			   int type3, int type4,
			   double target, double force) {
	this.type1 = type1;
	this.type2 = type2;
	this.type3 = type3;
	this.type4 = type4;
    	this.force = force;
	force2 = 2*force;
	this.target = Angle.deg2rad(target);
    }
 
    public int compareTo(Object obj) {
	OutOfPlaneParameters other = (OutOfPlaneParameters) obj;
	if (type1 > other.type1) return 1;
	if (type1 < other.type1) return -1;
	if (type2 > other.type2) return 1;
	if (type2 < other.type2) return -1;
	if (type3 > other.type3) return 1;
	if (type3 < other.type3) return -1;
	if (type4 > other.type4) return 1;
	if (type4 < other.type4) return -1;
	return 0;
    }
	
    public String toString() {
	return "OutOfPlaneParameters\n"+
	    "\t type1          = "+Atom.type(type1)+
	    "\t type2          = "+Atom.type(type2)+
	    "\t type3          = "+Atom.type(type3)+
	    "\t type4          = "+Atom.type(type4)+
	    "\t target   = "+Angle.rad2deg(target)+"\n"+
	    "\t force = "+force;
    }

    public Parameters create(StringTokenizer line) {
	return new OutOfPlaneParameters(line);
    }

    public Filter isA() {
	return new isA();
    }

    private  class isA implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof OutOfPlaneParameters);
	}
    }
}
