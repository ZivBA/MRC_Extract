package meshi.energy.caAttraction;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/** 
*
*  
**/

public class CaAttractionCreator extends EnergyCreator  implements KeyWords {
	
	/** 
	 **/ 
    public CaAttractionCreator(double weight) {
  		super(weight);
    }
	
    /** 
     * default constructor
     **/
    public CaAttractionCreator() {
  		super(1.0);
    }
	
	/** 
	 **/	    
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					                        CommandList commands) {
		return new CaAttractionEnergy(protein.atoms().CAFilter(), distanceMatrix, weight());
	}
    	
}//end
