package meshi.energy.bond;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.molecularElements.Atom;
import meshi.util.filters.Filter;
/**
 * The energy parameters associated with a bond between two atoms.
 * These parameters depend on the atom types.
 **/
public class BondParameters extends Parameters implements Comparable{
    /**
     * The target distance for this type of bond.
     **/ 
    public final double target;
    
    /**
     * The force constant for this type of bond.
     **/ 
    public final double force;
 
    /**
     * The square of the force constant.
     * Not a real parameter but saves time.
     **/ 
    public final double force2;

    /**
     * Atom type of the first atom.
     * The convention is that it is the smaller atom type in the pair.
     **/
    public final int type1;

    /**
     * Atom type of the second atom.
     * The convention is that it is the larger atom type in the pair.
     **/
     public final int type2;

    /**
     * Auxiliary constructor.
     **/
    public BondParameters() {this(-999,-999,-999,999);}

    /** 
     * Constructs BondParmeters from a string  (typically a line in the bondParameters file).
     **/
    public BondParameters(String line) {
	this(new StringTokenizer(line));
    }

    /**
     * A helper constructor for ( BondParameters(String line) ).
     **/
    public BondParameters(StringTokenizer line) {
	this(Atom.type(line.nextToken()), // type1
	     Atom.type(line.nextToken()), // type2
	     toDouble(line.nextToken()), // targetDistance
	     toDouble(line.nextToken())); // forceConstant
    }

    /**
     * A helper constructor for ( BondParameters(String line) ).
     **/
    public BondParameters(int type1, int type2, 
			  double  targetDistance, double forceConstant) {
	if (type1 < type2) {
	    this.type1 = type1;
	    this.type2 = type2;
	}
	else {
	    this.type1 = type2;
	    this.type2 = type1;
	}
	target = targetDistance ;
	force = forceConstant;
	force2 = 2.0 * forceConstant;
    }

    /**
     * For search keys in parameterList
     **/
    public BondParameters(int type1, int type2) {
	this(type1, type2, -1, -1);
    }
    

    /**
     * Defines order within objects and thus allows sorting for more efficient searches.
     **/
    public int compareTo(Object other) {
	if (! (other instanceof BondParameters))
	    throw new RuntimeException("Weird argument to "+
				       "BondParameters.compairTo(Object other)");
	BondParameters bp = (BondParameters) other;
	if (type1 > bp.type1) return 1;
	if (type1 < bp.type1) return -1;
	if (type2 > bp.type2) return 1;
	if (type2 < bp.type2) return -1;
	return 0;
    }
	
    public String toString() {
	return "BondParameters\n"+
	    "\t type1  = "+Atom.type(type1)+"\n"+
	    "\t type2  = "+Atom.type(type2)+"\n"+
	    "\t target = "+target+"\n"+
	    "\t force  = " +force;
    }

    private  class isA implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof BondParameters);
	}
    }
}
