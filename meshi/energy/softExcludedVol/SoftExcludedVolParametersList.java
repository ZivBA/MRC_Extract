package meshi.energy.softExcludedVol;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.geometry.Distance;
import meshi.molecularElements.Atom;
import meshi.util.file.MeshiLineReader;


public class SoftExcludedVolParametersList extends ParametersList {
	private double frac;

    public SoftExcludedVolParametersList(String parametersFileName) {this(parametersFileName,1.0);}

    	
    public SoftExcludedVolParametersList(String parametersFileName , double frac) {
    	super();
    	this.frac = frac;
    	SoftExcludedVolParameters tmpParam;
    	int maxType = -1;
    	SoftExcludedVolParameters[] tmpAr;
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
	    		tmpParam = (SoftExcludedVolParameters) createParameters(line);
	    		if (tmpParam.largeType>maxType)
	    		   maxType = tmpParam.largeType;
	    	}
	    }
	    catch (RuntimeException e) {
	    	System.out.println("A problem while reading parameters file "+
			   parametersFileName);
			   throw e;
		}
		tmpAr = new SoftExcludedVolParameters[((maxType+1)*(maxType+2))/2];
	    // Second reading actually creates the SoftExcludedVolParameters objects.
	    try {
	    	lines = new MeshiLineReader(parametersFileName);
	    	while ((line = lines.readLine("#")) != null) {
	    		tmpParam = (SoftExcludedVolParameters) createParameters(line);
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
	return new SoftExcludedVolParameters(new StringTokenizer(line) , frac);
    }
    
    public Parameters parameters(Object obj) {
	Distance distance = (Distance) obj;
	int largeType = distance.largeType();
	int smallType = distance.smallType();
	try {
	    return (SoftExcludedVolParameters) elementAt(largeType*(largeType+1)/2+smallType);
	}
	catch (Exception e) {
	    throw new RuntimeException(largeType+"-"+Atom.type(largeType)+" "+ 
				       smallType+"-"+Atom.type(smallType)+"\n"+e); 
	}	
   }
}
