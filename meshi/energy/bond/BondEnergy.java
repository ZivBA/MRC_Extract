package meshi.energy.bond;
import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.AtomPair;
import meshi.molecularElements.AtomPairList;
/**
 * Bond energy term. 
 * Has the general form <b> Eb = SIGMAi(Ki(Di-D0i)^2)</b>
 * where <b>Di</b> 
 * is the distance between the two bonded atoms, <b> D0i </b> 
 * is their expected average distance (depends on their types) and <b>Ki</b>
 * is a force constant that again, depends on the atom types.<br>
 * This class is used for both calculating the bond-energy term of an energy function 
 * and for updating the forces on each atom accordingly.<b>
 * It is assumed that the list of bonds is constant during the simulation. That is 
 * no bonds are made or broken.
 */
public class BondEnergy extends SimpleEnergyTerm {
    /**
     * The constructor associates any bond with its parameters.
     **/
    protected DistanceMatrix distanceMatrix;

    public BondEnergy() {}

    public BondEnergy(AtomPairList bondList, 
		      DistanceMatrix distanceMatrix,
		      BondParametersList  parametersList, 
		      double weight) {
	super(toArray(distanceMatrix), parametersList, weight);
	comment = "Bond";
	this.distanceMatrix = distanceMatrix;
	createElementsList(bondList);
    }

 
    public EnergyElement createElement(Object baseElement, Parameters parameters) {
	return new BondEnergyElement(((AtomPair)baseElement), parameters, distanceMatrix, weight);
    }

}
	
	
