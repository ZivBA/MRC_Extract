
package meshi.energy.LennardJones;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.geometry.Distance;
import meshi.molecularElements.Atom;
import meshi.util.file.MeshiLineReader;



public class LennardJonesParametersList extends ParametersList {
	private double frac;

    public LennardJonesParametersList(String parametersFileName) {this(parametersFileName,1.0);}

    	
    public LennardJonesParametersList(String parametersFileName , double frac) {
    	super();
    	this.frac = frac;
    	LennardJonesParameters tmpParam;
    	int maxType = -1;
    	LennardJonesParameters[] tmpAr;
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
	    		tmpParam = (LennardJonesParameters) createParameters(line);
	    		if (tmpParam.largeType>maxType)
	    		   maxType = tmpParam.largeType;
	    	}
	    }
	    catch (RuntimeException e) {
	    	System.out.println("A problem while reading parameters file "+
			   parametersFileName);
			   throw e;
		}
		tmpAr = new LennardJonesParameters[((maxType+1)*(maxType+2))/2];
	    // Second reading actually creates the LennardJonesParameters objects.
	    try {
	    	lines = new MeshiLineReader(parametersFileName);
	    	while ((line = lines.readLine("#")) != null) {
	    		tmpParam = (LennardJonesParameters) createParameters(line);
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
	return new LennardJonesParameters(new StringTokenizer(line) , frac);
    }
    
    public Parameters parameters(Object obj) {
	Distance distance = (Distance) obj;
	int largeType = distance.largeType();
	int smallType = distance.smallType();
	try {
	    return (LennardJonesParameters) elementAt(largeType*(largeType+1)/2+smallType);
	}
	catch (Exception e) {
	    throw new RuntimeException(largeType+"-"+Atom.type(largeType)+" "+ 
				       smallType+"-"+Atom.type(smallType)+"\n"+e); 
	}	
   }
}
