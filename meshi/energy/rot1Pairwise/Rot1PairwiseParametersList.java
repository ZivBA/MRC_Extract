
package meshi.energy.rot1Pairwise;

import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.geometry.Distance;
import meshi.molecularElements.Atom;
import meshi.util.file.MeshiLineReader;



public class Rot1PairwiseParametersList extends ParametersList {
    	
    public Rot1PairwiseParametersList(String parametersFileName) {
    	super();
    	Rot1PairwiseParameters tmpParam;
    	Rot1PairwiseParameters[] tmpAr;
    	MeshiLineReader lines;
    	int maxType = -1;
    	String line;
    	
	    if (parametersFileName == null) 
	    	throw new RuntimeException("No parameter file name in " + this);
	    	 
	    System.out.println("Loading "+this+" parameters from "+parametersFileName);
	    // First reading to see that everything is ok, and also to see what is the 
	    // largest atom type.
	    try {
	    	lines = new MeshiLineReader(parametersFileName);
	    	while ((line = lines.readLine("#")) != null) {
	    		tmpParam = (Rot1PairwiseParameters) createParameters(line);
	    		if (tmpParam.largerType>maxType)
	    		   maxType = tmpParam.largerType;
	    	}
	    }
	    catch (RuntimeException e) {
	    	System.out.println("A problem while reading parameters file "+
			   parametersFileName);
			   throw e;
		}
		tmpAr = new Rot1PairwiseParameters[((maxType+1)*(maxType+2))/2];
	    // Second reading actually creates the Rot1PairwiseParameters objects.
	    try {
	    	lines = new MeshiLineReader(parametersFileName);
	    	while ((line = lines.readLine("#")) != null) {
	    		tmpParam = (Rot1PairwiseParameters) createParameters(line);
	    		tmpAr[(tmpParam.largerType*(tmpParam.largerType+1))/2+tmpParam.smallerType] = tmpParam;
	    	}
	    }
	    catch (RuntimeException e) {
	    	System.out.println("A problem while reading parameters file "+
			   parametersFileName);
			   throw e;
		}
		
		for (int c=0 ; c<tmpAr.length ; c++)
		   add(tmpAr[c]);
	}

    public Parameters createParameters(String line) {
    	return new Rot1PairwiseParameters(new StringTokenizer(line));
    }
    
    public Parameters parameters(Object obj) {
    	Distance distance = (Distance) obj;
    	int largeType = distance.largeType();
    	int smallType = distance.smallType();
    	try {
    		return (Rot1PairwiseParameters) elementAt(largeType*(largeType+1)/2+smallType);
    	}
    	catch (Exception e) {
    		throw new RuntimeException(largeType+"-"+Atom.type(largeType)+" "+ 
    				smallType+"-"+Atom.type(smallType)+"\n"+e); 
    	}	
    }
}
