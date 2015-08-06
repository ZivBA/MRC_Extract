package meshi.energy.rot1Pairwise;

import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.molecularElements.Atom;
import meshi.util.filters.Filter;

public class Rot1PairwiseParameters extends Parameters{
	public final int smallerType; 
	public final int largerType;
	public final double[] EofR;
	public double maxDis;

    public Rot1PairwiseParameters() {
    	EofR = null;
    	smallerType = largerType = -1;
    	maxDis = -1.0;
    }

    public Rot1PairwiseParameters(StringTokenizer st) {
    	int first = Atom.type(st.nextToken());
    	int second = Atom.type(st.nextToken());
    	if (first>second) {
    	    smallerType = second;
    	    largerType = first;
    	}
    	else {
    	    smallerType = first;
    	    largerType = second;
    	}	
    	EofR = new double[st.countTokens()];
    	for (int counter=0; st.hasMoreTokens() ; counter++) {
    		EofR[counter] = (new Double(st.nextToken())).doubleValue();
    	}
    	maxDis = EofR.length-1;
    }
    
    public Filter isA() {return (new isA());}

    public Parameters create(StringTokenizer stringTokenizer) {
    	return (new  Rot1PairwiseParameters(stringTokenizer));
    }

    private  class isA implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof Rot1PairwiseParameters);
	}
    }

    public String toString() {
    	return smallerType+"\t"+largerType;
    }
}


