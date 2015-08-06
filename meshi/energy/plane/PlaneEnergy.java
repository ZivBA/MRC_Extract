package meshi.energy.plane;
import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionList;
/**
 * Plane energy term. 
 **/
public class PlaneEnergy extends SimpleEnergyTerm{
    protected TorsionList torsionList;
    protected DistanceMatrix distanceMatrix;

    public PlaneEnergy() {}

    public PlaneEnergy(TorsionList torsionList, DistanceMatrix distanceMatrix, 
		       PlaneParametersList  parametersList, double weight) {
	super(toArray(distanceMatrix, torsionList), parametersList, weight);
	this.torsionList = torsionList;
	this.distanceMatrix = distanceMatrix;
	createElementsList(torsionList);
	comment = "Plane";
    }




    public EnergyElement createElement(Object baseElement, Parameters parameters) {
	return new PlaneEnergyElement(((Torsion)baseElement), parameters, weight);
    }
}    
