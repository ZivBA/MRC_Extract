package meshi.energy.torsionVal;
import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionList;
/**
 * TorsionVal energy term. 
 **/
public class TorsionValEnergy extends SimpleEnergyTerm{

    public TorsionValEnergy() {}

    public TorsionValEnergy(DistanceMatrix distanceMatrix,
			    TorsionList torsionList, 
			    TorsionValParametersList  parametersList, 
			    double weight) {
	super(toArray(distanceMatrix, torsionList), parametersList, weight);
	createElementsList(torsionList);
	comment = "TorsionVal";
    }
 

    public EnergyElement createElement(Object baseElement, Parameters parameters) {
	return new TorsionValEnergyElement(((Torsion)baseElement), parameters, weight);
    }


    public void handleMissingParameters(Object obj) {}
    
    public TorsionValEnergyElement getTorsionEnergyElement(int residueNumber, String name) {
    	for (int c=0 ; c<elementsList.size() ; c++)
    		if ((((TorsionValEnergyElement) elementsList.elementAt(c)).torsion().getTorsionResNum()==residueNumber) &&
    				((TorsionValEnergyElement) elementsList.elementAt(c)).torsion().getTorsionName().equals(name))
    			return ((TorsionValEnergyElement) elementsList.elementAt(c));
    	return null;
    }
    
    public void setChisTargetsAutomaticly() {
    	for (int ind=0 ; ind<elementsList().size() ; ind++) {
    		TorsionValEnergyElement torElement = (TorsionValEnergyElement) elementsList().elementAt(ind);
    		torElement.setTarget(torElement.torsion().torsion());
    	}
    }

}

	
