package meshi.energy.distanceConstrains;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.util.filters.Filter;
/**
 * The energy parameters associated with a distanceConstrains between two atoms.
 **/
public class DistanceConstrainParameters extends Parameters implements Comparable{
    /**
     * The target distance for this distanceConstrain.
     **/ 
    public final double target;
    
    /**
     * The tolerance this distanceConstrain.
     **/ 
    public final double tolerance;
    
    /**
     * The force constant for this type of distanceConstrain.
     **/ 
    public final double force;
 
    /**
     * The square of the force constant.
     * Not a real parameter but saves time.
     **/ 
    public final double force2;

    /**
     * Residue name of the first atom.
     **/
    public final String residue1Name;

    /**
     * Residue number of the first atom.
     **/
    public final int residue1Number;

    /**
     * Name of the the first atom.
     **/
    public final String name1;

    /**
     * Residue name of the second atom.
     **/
    public final String residue2Name;

    /**
     * Residue number of the second atom.
     **/
    public final int residue2Number;

    /**
     * Name of the the second atom.
     **/
    public final String name2;


    public DistanceConstrainParameters(String line) {
	this(new StringTokenizer(line));
    }
    
    public DistanceConstrainParameters(StringTokenizer line) {
	this(line.nextToken(),             // residue1 name
	     toInt(line.nextToken()),      // esidue1 number
	     line.nextToken(),             // atom1 name
	     line.nextToken(),             // residue2 name
	     toInt(line.nextToken()),      // esidue2 number
	     line.nextToken(),             // atom2 name
	     toDouble(line.nextToken()),   // targetDistance
	     toDouble(line.nextToken()),   // tolerance
	     toDouble(line.nextToken()));  // forceConstant
    }

    public DistanceConstrainParameters(String residue1Name, int residue1Number, String name1,
				       String residue2Name, int residue2Number, String name2,
				       double target, double tolerance, double force) {
	
	this.residue1Name = residue1Name;
	this.residue1Number = residue1Number;
	this.name1 = name1;
	this.residue2Name = residue2Name;
	this.residue2Number = residue2Number;
	this.name2 = name2;

	this.target = target;
	this.tolerance = tolerance;
	this.force = force;
	force2 = 2.0 * force;
    }

    private  class isA implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof DistanceConstrainParameters);
	}
    }
    public int compareTo(Object obj) {
	throw new RuntimeException("compareTo not applicable");
    }

    public String toString() {
	return ("DistanceConstrainParameters: "+residue1Name+"_"+residue1Number+"_"+name1+
		" "+residue2Name+"_"+residue2Number+"_"+name2+" "+target+" "+tolerance+" "+force);
    }

    public String residue1Name() {
	return residue1Name;
    }
    public String residue2Name() {
	return residue2Name;
    }
    public String name1() {
	return name1;
    }
    public String name2() {
	return name2;
    }
    public int residue1Number() {
	return residue1Number;
    }
    public int residue2Number() {
	return residue2Number;
    }
}
