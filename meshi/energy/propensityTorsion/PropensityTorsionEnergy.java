package meshi.energy.propensityTorsion;
import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.TorsionPair;
import meshi.geometry.TorsionPairList;
import meshi.util.filters.Filter;

/**
 *This class minimizes the propensity of the PHI,PSI residue over the entire protein. It works practically the same 
 *as TwoTorsionsParameters.java with minor variations:
 *1) PropensityTorsionEnergyElement.java calculates propensity instead of probability.
 *2) The secondary structure is disregarded so that the parameter file format is slightly changed (the SS field is omited)
 **/
 
public class PropensityTorsionEnergy extends SimpleEnergyTerm{
    public static final Filter isTorsionPair = new IsTorsionPair();
    protected TorsionPairList torsionPairList;
    protected DistanceMatrix distanceMatrix;

    public PropensityTorsionEnergy() {}

    public PropensityTorsionEnergy(TorsionPairList torsionPairList, DistanceMatrix distanceMatrix, 
		       PropensityTorsionParametersList  parametersList, double weight) {
	super(toArray(distanceMatrix, torsionPairList), parametersList, weight);
	this.torsionPairList = torsionPairList;
	this.distanceMatrix = distanceMatrix;
	createElementsList(torsionPairList);
	comment = "propensityTorsion";
    }




    public EnergyElement createElement(Object baseElement, Parameters parameters) {
	return new PropensityTorsionEnergyElement((TorsionPair) baseElement,
	                                    (PropensityTorsionParameters) parameters, weight);
    }

    public void handleMissingParameters(Object obj) {}


    public static class IsTorsionPair implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof TorsionPair);
	}
    }


}    
