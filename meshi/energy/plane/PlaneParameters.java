package meshi.energy.plane;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.molecularElements.Atom;
import meshi.util.filters.Filter;
public class PlaneParameters extends Parameters implements Comparable {
    public final static int TRANS = 1;
    public final static int CIS = 2;
    public final static int CIS_TRANS = 3;

    public final double force, force2;
    public int trans;
    public final int type1, type2, type3, type4;
    static String temp;
    public PlaneParameters() {
	this(-1,-1,-1,-1,-9999.9,1); 
    }
    public PlaneParameters(int type1, int type2, int type3, int type4) {
	this(type1,type2,type3,type4,-9999.9,1); 
    }
    public PlaneParameters(String line) {
	this(new StringTokenizer(line));
    }
    private PlaneParameters(StringTokenizer line) {
	this(Atom.type(line.nextToken()), 
	     Atom.type(line.nextToken()), 
	     Atom.type(line.nextToken()), 
	     Atom.type(line.nextToken()), 
	     toDouble(line.nextToken()), 
	     (temp = line.nextToken()).equals("TRANS")?TRANS:
	                                      (temp.equals("CIS")?CIS:CIS_TRANS));
    }

    public PlaneParameters(int type1, int type2,
			   int type3, int type4,
			   double force, int trans) {
	this.type1 = type1;
	this.type2 = type2;
	this.type3 = type3;
	this.type4 = type4;
	this.force = force;
	force2 = 2*force;
	this.trans = trans;
    }
 
    public int compareTo(Object obj) {
	PlaneParameters other = (PlaneParameters) obj;
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
	return "PlaneParameters\n"+
	    "\t type1   = "+Atom.type(type1)+
	    "\t type2   = "+Atom.type(type2)+
	    "\t type3   = "+Atom.type(type3)+
	    "\t type4   = "+Atom.type(type4)+
	    "\t trans   = "+trans+"\n"+
	    "\t force = "+force;
    }

    public Parameters create(StringTokenizer line) {
	return new PlaneParameters(line);
    }

    public Filter isA() {
	return new isA();
    }

    private  class isA implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof PlaneParameters);
	}
    }
}
