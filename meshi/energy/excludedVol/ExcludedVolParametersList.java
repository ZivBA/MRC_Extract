package meshi.energy.excludedVol;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.geometry.Distance;
import meshi.molecularElements.Atom;
import meshi.util.file.MeshiLineReader;

/**
* The elements in the list of ExcludedVolParametersList are sorted by the atom types in the 
* following way:
* [EV parameters of type 0 with type 0]
* [EV parameters of type 1 with type 0]
* [EV parameters of type 1 with type 1]
* [EV parameters of type 2 with type 0]
* [EV parameters of type 2 with type 1]
* [EV parameters of type 2 with type 2]
* [EV parameters of type 3 with type 0]
* [EV parameters of type 3 with type 1]
* [EV parameters of type 3 with type 2]
* etc.
*
**/

public class ExcludedVolParametersList extends ParametersList {
	
    public ExcludedVolParametersList(String parametersFileName) {
    	super();
    	
    	ExcludedVolParameters tmpParam;
    	int maxType = -1;
    	ExcludedVolParameters[] tmpAr;
    	MeshiLineReader lines;
    	String line;
    	
	    if (parametersFileName == null) 
	    	throw new RuntimeException("No parameter file name in " + this);
	    	 
	    System.out.println("Loading "+this+" parameters from "+parametersFileName);
	    // First reading to see that everything is ok, and also to see what is the 
	    // largest atom type.
	    try {
	    	lines = new MeshiLineReader(parametersFileName);
	    	while ((line = lines.readLine("#")) != null) {
	    		tmpParam = (ExcludedVolParameters) createParameters(line);
	    		if (tmpParam.largeType>maxType)
	    		   maxType = tmpParam.largeType;
	    	}
	    }
	    catch (RuntimeException e) {
	    	System.out.println("A problem while reading parameters file "+
			   parametersFileName);
			   throw e;
		}
		tmpAr = new ExcludedVolParameters[((maxType+1)*(maxType+2))/2];
	    // Second reading actually creates the ExcludedVolParameters objects.
	    try {
	    	lines = new MeshiLineReader(parametersFileName);
	    	while ((line = lines.readLine("#")) != null) {
	    		tmpParam = (ExcludedVolParameters) createParameters(line);
	    		tmpAr[tmpParam.largeType*(tmpParam.largeType+1)/2+tmpParam.smallType] = tmpParam;
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
	return new ExcludedVolParameters(new StringTokenizer(line));
    }
    
    public Parameters parameters(Object obj) {
	Distance distance = (Distance) obj;
	int largeType = distance.largeType();
	int smallType = distance.smallType();
	try {
	    return (ExcludedVolParameters) elementAt(largeType*(largeType+1)/2+smallType);
	}
	catch (Exception e) {
	    throw new RuntimeException(largeType+"-"+Atom.type(largeType)+" "+ 
				       smallType+"-"+Atom.type(smallType)+"\n"+e); 
	}	
   }
}
