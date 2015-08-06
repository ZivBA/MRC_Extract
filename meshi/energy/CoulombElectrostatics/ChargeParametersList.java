package meshi.energy.CoulombElectrostatics;

import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.molecularElements.Atom;
import meshi.util.file.MeshiLineReader;

/**
* 
*
* A list that holds all the electrical charges of all atoms in all the amino acids
* (ChargeParameters).
* The charge data is taken from the data file: "CHARGE_PARAMETERS".
*
* Constructors
*	public ChargeParametersList(String parametersFileName)	-  takes the data file name as a parameter.
*	
* Object methods 
*	public Parameters createParameters(String line)
*   public Parameters getParameters(Object obj)
*	public String toString()
*
*  
**/
public class ChargeParametersList extends ParametersList {
	private final ChargeParameter zero = new ChargeParameter(-1,0);
	/**
	 * A constructor that takes the charge's data file name. It creates a
	 * list that holds objects of the type "ChargeParameter".
	 * @param parametersFileName
	 **/
	public ChargeParametersList(String parametersFileName) {
	    super();
		ChargeParameter[] tmpAr;
    	MeshiLineReader lines;
    	String line;
    	ChargeParameter tmpParam;
    	int maxType = -1;
    	if (parametersFileName == null) 
   	    	throw new RuntimeException("No parameter file name in " + this);
   	    	 
   	    System.out.println("Loading "+this+" parameters from "+parametersFileName);
   	    try {
	    	lines = new MeshiLineReader(parametersFileName);
	    	while ((line = lines.readLine("#")) != null) {
	    		tmpParam = (ChargeParameter) createParameters(line);
	    		if (tmpParam.atomType>maxType)
	    		   maxType = tmpParam.atomType;
	    	}
	    }
	    catch (RuntimeException e) {
	    	System.out.println("A problem while reading parameters file "+
			   parametersFileName);
			   throw e;
		}
		tmpAr = new ChargeParameter[maxType+1];
		this.size=tmpAr.length;
	    // Second reading actually creates the ExcludedVolParameters objects.
	    try {
	    	lines = new MeshiLineReader(parametersFileName);
	    	while ((line = lines.readLine("#")) != null) {
	    		tmpParam = (ChargeParameter) createParameters(line);
	    		tmpAr[tmpParam.atomType] = tmpParam;
	    	}
	    }
	    catch (RuntimeException e) {
	    	System.out.println("A problem while reading parameters file "+
			   parametersFileName);
			   throw e;
		}
		
		for (int c=0 ; c<tmpAr.length ; c++){
		   if (tmpAr[c]!=null) set(tmpAr[c].atomType,tmpAr[c]);		   
		}
	}

   
	

	/**
	 * Each line in the file represents a ChargeParameter of an atom. This function takes a line and 
	 * calls the Constructor of ChargeParameter, which in turn builds a new Parameter object form that line.
	 * @param line
	 * @return Parameters
	 **/
	
	public Parameters createParameters(String line) {
		return new ChargeParameter(new StringTokenizer(line));
    }

	/**
	 * This function takes an Atom object and returns its ChargeParameter.
	 * i.e. the relevant charge for the atom in its position in the amino acid.
	 * @param obj a specific Atom
	 * @return Parameters
	 **/ 
    public Parameters parameters(Object obj) {
    	Atom atom = (Atom) obj;
   		int atomType = atom.type;
   		
		try {
		//	ChargeParameter  charge_parameter=new ChargeParameter(atomType,0);			
			ChargeParameter  charge_parameter;
			if (this.internalArray[atomType]!=null)			
		    charge_parameter = (ChargeParameter) elementAt(atomType);
		//	else charge_parameter=new ChargeParameter(atomType,0);
			else charge_parameter=zero;
			return charge_parameter;
		}
		catch (Exception e) {
			
		    throw new RuntimeException("atomType " +"- "+ atomType + "\n" + e); 
		}	
   }
    
    /**
     * toString
     * @return a String representation of this object
     **/	
	public String toString(){
     	String ans = "" ;
        Object []  obj = this.toArray();
        for(int i =0; i < obj.length; i++){
			ChargeParameter  cp = (ChargeParameter) obj [i];
			ans = ans + i + " " + cp.toString();
			ans = ans + "\n";
		}
		return ans;
	}
}
