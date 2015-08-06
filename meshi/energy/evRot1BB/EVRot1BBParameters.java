package meshi.energy.evRot1BB;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.molecularElements.Atom;
import meshi.util.filters.Filter;

public class EVRot1BBParameters extends Parameters {
    public final double  sigma;
    public final int smallType;
    public final int largeType;
    
    
    public EVRot1BBParameters() {
    smallType = -1;
    largeType = -1;	
	sigma = -1;
    }

    public EVRot1BBParameters(StringTokenizer st) {
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
    }
    
    public Filter isA() {return (new isA());}

    public Parameters create(StringTokenizer stringTokenizer) {
	return (new  EVRot1BBParameters(stringTokenizer));
    }

    private  class isA implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof EVRot1BBParameters);
	}
    }

    public String toString() {return "largeType:" + Atom.type(largeType)+"   "+ 
				       "smallType:"+Atom.type(smallType)+"   sigma:"+sigma + "\n";}
}
