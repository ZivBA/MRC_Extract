package meshi.energy.disConst;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.util.MeshiList;


public class DisConstEnergy  extends SimpleEnergyTerm{
    /**
     * The constructor associates any bond with its parameters.
     **/

    public DisConstEnergy() {}

    public DisConstEnergy(DistanceList disList, double weight) {
	
	elementsList = new MeshiList();
	this.weight = weight;

	for (int disCounter=0 ; disCounter<disList.size(); disCounter++) {
		Distance dis = disList.distanceAt(disCounter);
	    EnergyElement newElement = createElement(dis);
	    if (!newElement.frozen())
	    	elementsList.add(newElement);	    
    }
	comment = "DisConst";
    }


    public EnergyElement createElement(Object baseElement) {
    	return new DisConstEnergyElement(weight, (Distance) baseElement);
    }

     public void update(){}
     public void update(int i ){}

     public EnergyElement createElement(Object baseElement,Parameters p){return null;}

}


