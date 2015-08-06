package meshi.energy.LennardJones;

import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.molecularElements.Atom;
import meshi.util.filters.Filter;

public class LennardJonesParameters extends Parameters{
    public final double epsilon = 0.15; // Average of values for C,N and O in AMBER94
	public final double  sigma;
	public final int smallType;
	public final int largeType;
    public final double sigma6, sigma6EpsilonFour, minusTwelveSigma6;  

    public LennardJonesParameters() {
    	sigma = sigma6 = sigma6EpsilonFour = minusTwelveSigma6 = -1;
    	smallType = largeType = -1;
    }

    public LennardJonesParameters(StringTokenizer st) {
    	this(st,1.0);
    }

    public LennardJonesParameters(StringTokenizer st,double frac) {
    	int first = Atom.type(st.nextToken());
    	int second = Atom.type(st.nextToken());
    	if (first>second) {
    	    smallType = second;
    	    largeType = first;
    	}
    	else {
    	    smallType = first;
    	    largeType = second;
    	}	
    	sigma = toDouble(st.nextToken())*frac;
    	sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
    	sigma6EpsilonFour = 4*epsilon*sigma6;
    	minusTwelveSigma6 = -12*sigma6;
    }
    
    public Filter isA() {return (new isA());}

    public Parameters create(StringTokenizer stringTokenizer) {
    	return (new  LennardJonesParameters(stringTokenizer));
    }

    private  class isA implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof LennardJonesParameters);
	}
    }

    public String toString() {
    	return epsilon+"\t"+sigma;
    }
}


