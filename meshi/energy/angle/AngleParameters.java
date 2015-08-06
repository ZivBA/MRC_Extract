package meshi.energy.angle;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.geometry.Angle;
import meshi.molecularElements.Atom;
import meshi.util.filters.Filter;
public class AngleParameters extends Parameters implements Comparable {
    public final double target, force, force2;
    public final int type1, type2, type3;
    public AngleParameters() {
	this(-1,-1,-1,-9999.9,-9999.9); 
    }
    public AngleParameters(int type1, int type2, int type3) {
	this(type1,type2,type3,-9999.9,-9999.9); 
    }
    public AngleParameters(String line) {
	this(new StringTokenizer(line));
    }
    private AngleParameters(StringTokenizer line) {
	this(Atom.type(line.nextToken()), 
	     Atom.type(line.nextToken()), 
	     Atom.type(line.nextToken()), 
	     toDouble(line.nextToken()), toDouble(line.nextToken()));
    }

    public AngleParameters(int type1, int type2, int type3, double  targetAngle, double forceConstant) {
	this.type2 = type2;
	if (type1 < type3) {
	    this.type1 = type1;
	    this.type3 = type3;
	}
	else {
	    this.type1 = type3;
	    this.type3 = type1;
	}
	target= Angle.deg2rad(targetAngle);
	force = forceConstant;
	force2 = 2.0 * force;
    }
    public int compareTo(Object other) {
	if (! (other instanceof AngleParameters))
	    throw new RuntimeException("Weird argument to "+
				       "AngleParameters.compairTo(Object other)");
	AngleParameters aep = (AngleParameters) other;
	if (type1 > aep.type1) return 1;
	if (type1 < aep.type1) return -1;
	if (type2 > aep.type2) return 1;
	if (type2 < aep.type2) return -1;
	if (type3 > aep.type3) return 1;
	if (type3 < aep.type3) return -1;
	return 0;
    }
	
    public String toString() {
	return "AngleParameters\n"+
	    "\t type1          = "+Atom.type(type1)+
	    "\t type2          = "+Atom.type(type2)+
	    "\t type3          = "+Atom.type(type3)+
	    "\t target   = "+Angle.rad2deg(target)+"\n"+
	    "\t force = "+force;
    }

    public Parameters create(StringTokenizer line) {
	return new AngleParameters(line);
    }

    public Filter isA() {
	return new isA();
    }

    private  class isA implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof AngleParameters);
	}
    }
}
