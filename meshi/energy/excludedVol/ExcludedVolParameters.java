package meshi.energy.excludedVol;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.molecularElements.Atom;
import meshi.util.filters.Filter;

public class ExcludedVolParameters extends Parameters {
    public final double  sigma;
    public final int smallType;
    public final int largeType;
    /** 
    * ALPHA is the transition zone (in Ang) where the energy change in the forth power 0.0 to 1.0.
    **/    
    public final double  ALPHA = 0.2; 
    /**
    * C is the parameter of  EV = C*(d-sigma)^4 in the range [0,sigma]
    **/    
    public final double  C;
    
    
    public ExcludedVolParameters() {
    smallType = -1;
    largeType = -1;	
	sigma = -1;
	C =  -1;
    }

    public ExcludedVolParameters(StringTokenizer st) {
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
	sigma = toDouble(st.nextToken());
	C =  1/(ALPHA*ALPHA*ALPHA*ALPHA);
    }
    
    public Filter isA() {return (new isA());}

    public Parameters create(StringTokenizer stringTokenizer) {
	return (new  ExcludedVolParameters(stringTokenizer));
    }

    private  class isA implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof ExcludedVolParameters);
	}
    }

    public String toString() {return ""+sigma;}
}
