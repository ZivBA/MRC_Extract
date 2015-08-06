package meshi.energy.CoulombElectrostatics;

import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.molecularElements.Atom;

/**
* 
*
* Electric charge of an atom.
* Every atom in every Amino acid has it's own charge.
* The data of the charges is taken from the data file: "CHARGE_PARAMETERS ".
*
* Object variables   
*	public final double  charge                     - the charge of the atom.
*
* Constructors
*	public  ChargeParameter()						- Builds a default ChargeParameter
*	public ChargeParameter(StringTokenizer st)		- Builds a ChargeParameter that takes a StringTokaenizer that holds 
*								 	                  a line from the data file, and extracts the charge from it.
*   public ChargeParameter(double charge)           - Builds a ChargeParameter that takes the charge value as a  parameter.
*
* Object methods 
*	public double charge()
*	public String toString()
*
*  
**/

public class ChargeParameter extends Parameters {
    public double charge; //charge   
    public int atomType;
    
    /**
     *  default ChargeParameter constructor.
     **/
    public ChargeParameter() {
		 charge = -1;
    }
    
    /**
     * A Constructor of ChargeParameter that takes a StringTokaenizer as its parameter.
     * The StringTokenizer holds a line from the data file and the charge is extracted from it.
     * @param st holds the relevant line from the data file
     **/ 
    public ChargeParameter(StringTokenizer st) {
	    String temp = st.nextToken();	    	  
	    atomType=Atom.type(temp);
	    temp=st.nextToken();	    
		charge = toDouble(temp);
    }
    /**
     * A Constructor of ChargeParameters that takes a charge value as a  parameter.
     * @param charge the charge value
     **/
    public ChargeParameter(int type,double charge){
    	atomType=type;
    	this.charge = charge;
    }
    
    /**
     * Returns the charge field
     * @returns the charge value of the chargeParameter
     **/
    public double charge(){
    	return charge;
    }
    
    /**
     * toString Function.
     * @return a String representation of this object
     **/	
    public String toString() {
    	return "type:\t"+atomType+"charge:\t"+ charge ;
    }

}//end

