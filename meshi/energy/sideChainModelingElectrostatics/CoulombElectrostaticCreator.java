
package meshi.energy.sideChainModelingElectrostatics;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/** 
* 
*
* Creates a CoulombElectrostaticTerm for the protein with a given weight.
* The creator: 
*  1. Reads data from the file, and creates the charge parameters list.
*  2. Initializes the DielectricConstsnt from the command file (given by user)
*  3. Creates a CoulombElectrostatics object that will evaluate the total electrostatics energy.
*	
* Constructors
*	public CoulombElectrostaticCreator(double weight)	-  A CoulombElectrostaticCreator that takes a weight
*                                                          parameter for energy term.                                                       
*	public CoulombElectrostaticCreator()    		    -  The default CoulombElectrostaticCreator.
*
* Object methods 
*	AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, CommandList commands)
*   public double getDielectricConstsnt(CommandList commands)
*
*  
**/

public class CoulombElectrostaticCreator extends EnergyCreator  implements KeyWords {
	
	/** 
	 * A CoulombElectrostaticCreator constructor that takes the weight
	 * given to this energy term as a parameter. 
	 * @param weight the weight given to this energy term.
	 **/ 
    public CoulombElectrostaticCreator(double weight) {
  		super(weight);
    }
	
    /** 
     * default constructor
     **/
    public CoulombElectrostaticCreator() {
  		super(ELECTROSTATICS);
    }
	
	/** 
	 * the creator: 
	 * 1. Reads data from the file "CHARGE_PARAMETERS" and creates the charge parameters list.
	 * 2. Initialize the DielectricConstsnt from the command file (given by user).
	 * 3. Creates a CoulombElectrostatics object that will evaluate the total electrostatics energy.
	 * @param protein
	 * @param distanceMatrix
	 * @param commands
	 * @return AbstractEnergy
	 **/	    
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					                        CommandList commands) {
		if (parametersList== null)
			parametersList = new ChargeParametersList(parametersDirectory(commands)+
							    "/"+ELECTROSTATICS_PARAMETERS);
		double dielectricConstant = 1.0;
		return new CoulombElectrostatics(distanceMatrix,(ChargeParametersList) parametersList, weight(), dielectricConstant);
	}
	
}//end
