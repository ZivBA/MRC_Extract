package meshi.energy.evRot1BB;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.geometry.Distance;
import meshi.molecularElements.Atom;
import meshi.util.file.MeshiLineReader;


public class EVRot1BBParameterList extends ParametersList {
	
	public final double maxEnd;
	
    public EVRot1BBParameterList(String parametersFileName) {
    	super();
    	
    	EVRot1BBParameters tmpParam;
    	int maxType = -1;
    	EVRot1BBParameters[] tmpAr;
    	MeshiLineReader lines;
    	String line;
    	double maxSig = -999999;
    	
	    if (parametersFileName == null) 
	    	throw new RuntimeException("No parameter file name in " + this);
	    	 
	    System.out.println("Loading "+this+" parameters from "+parametersFileName);
	    // First reading to see that everything is ok, and also to see what is the 
	    // largest atom type.
	    try {
	    	lines = new MeshiLineReader(parametersFileName);
	    	while ((line = lines.readLine("#")) != null) {
	    		tmpParam = (EVRot1BBParameters) createParameters(line);
	    		if (tmpParam.largeType>maxType)
	    		   maxType = tmpParam.largeType;
	    	}
	    }
	    catch (RuntimeException e) {
	    	System.out.println("A problem while reading parameters file "+
			   parametersFileName);
			   throw e;
		}
		tmpAr = new EVRot1BBParameters[((maxType+1)*(maxType+2))/2];
	    // Second reading actually creates the SoftExcludedVolParameters objects.
	    try {
	    	lines = new MeshiLineReader(parametersFileName);
	    	while ((line = lines.readLine("#")) != null) {
	    		tmpParam = (EVRot1BBParameters) createParameters(line);
	    		tmpAr[tmpParam.largeType*(tmpParam.largeType+1)/2+tmpParam.smallType] = tmpParam;
	    	}
	    }
	    catch (RuntimeException e) {
	    	System.out.println("A problem while reading parameters file "+
			   parametersFileName);
			   throw e;
		}
		
		for (int c=0 ; c<tmpAr.length ; c++) {
		   add(tmpAr[c]);
		   if (tmpAr[c].sigma>maxSig)
		   		maxSig = tmpAr[c].sigma;
		}
		maxEnd = maxSig;
	}

    public Parameters createParameters(String line) {
	return new EVRot1BBParameters(new StringTokenizer(line));
    }
    
    public Parameters parameters(Object obj) {
	Distance distance = (Distance) obj;
	int largeType = distance.largeType();
	int smallType = distance.smallType();
	try {
	    return (EVRot1BBParameters) elementAt(largeType*(largeType+1)/2+smallType);
	}
	catch (Exception e) {
	    throw new RuntimeException(largeType+"-"+Atom.type(largeType)+" "+ 
				       smallType+"-"+Atom.type(smallType)+"\n"+e); 
	}	
   }
}
