package meshi.energy.outOfPlane;
import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionList;
/**
 * OutOfPlane energy term. 
 **/
public class OutOfPlaneEnergy extends SimpleEnergyTerm{
    /**
     * The constructor associates any outOfPlane with its parameters.
     **/
    protected TorsionList torsionList;
    protected DistanceMatrix distanceMatrix;

    public OutOfPlaneEnergy() {}

    public OutOfPlaneEnergy(DistanceMatrix distanceMatrix,
			    TorsionList torsionList, 
			    OutOfPlaneParametersList  parametersList, 
			    double weight) {
	super(toArray(distanceMatrix, torsionList), parametersList, weight);
	this.torsionList = torsionList;
	this.distanceMatrix = distanceMatrix;
	createElementsList(torsionList);
	comment = "OutOfPlane";
    }


    public Parameters createParameters(String line) {
	return new OutOfPlaneParameters(line);
    }
 

    public EnergyElement createElement(Object baseElement, Parameters parameters) {
	return new OutOfPlaneEnergyElement(((Torsion)baseElement), parameters, weight);
    }

    public Parameters getKey(Object baseElement) {
	Torsion torsion = (Torsion) baseElement;
	return new OutOfPlaneParameters(torsion.atom1.type, torsion.atom2.type, 
					torsion.atom3.type, torsion.atom4.type);
    }

    public void handleMissingParameters(Object obj) {}
}	
	
