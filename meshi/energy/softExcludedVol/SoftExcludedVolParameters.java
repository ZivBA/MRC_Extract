package meshi.energy.softExcludedVol;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.molecularElements.Atom;
import meshi.util.filters.Filter;

public class SoftExcludedVolParameters extends Parameters {
    public final double  sigma;
    public final int smallType;
    public final int largeType;
    /** 
    * ALPHA is the transition zone (in Ang) where the energy change in the forth power 0.0 to 1.0.
    **/    
    public final double  ALPHA = 1.0; 
    /**
    * C is the parameter of  EV = C*(d-sigma)^4 in the range [0,sigma]
    **/    
    public final double  C;
    
    
    public SoftExcludedVolParameters() {
    smallType = -1;
    largeType = -1;	
	sigma = -1;
	C =  -1;
    }

    public SoftExcludedVolParameters(StringTokenizer st) {this(st,1.0);}

    public SoftExcludedVolParameters(StringTokenizer st,double frac) {
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
	C =  1/(ALPHA*ALPHA);
    }
    
    public Filter isA() {return (new isA());}

    public Parameters create(StringTokenizer stringTokenizer) {
	return (new  SoftExcludedVolParameters(stringTokenizer));
    }

    private  class isA implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof SoftExcludedVolParameters);
	}
    }

    public String toString() {return ""+sigma;}
}
