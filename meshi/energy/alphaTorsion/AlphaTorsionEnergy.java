package meshi.energy.alphaTorsion;
import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionList;

/**
 *This energy limits the alpha torsion (4 consecutive CAs) to be in a specific range.
 *It is also secondary structure sensitive, and the range depends on the residue SS.
 *It operates currently only on HELIX,SHEET secondary structure states, since the COIL,ALL 
 *states practically don't have any limitation on the alpha torsion.
 *
 *Important Note: This energy term must be accompanied by an ALPHA-angle energy term. This is 
 *because it has a non-continous point at torsion values of -Pi or Pi , and also when 
 *one of the 2 angles that make up the torsion is close to 0 or Pi. These discontinuites should 
 *not affect normal operation if the ALPHA-angle term is working. On very rare starting condition
 *these problems might never the less be encountered.  
 **/
 
             
public class AlphaTorsionEnergy extends SimpleEnergyTerm{
    protected TorsionList torsionList;
    protected DistanceMatrix distanceMatrix;

    public AlphaTorsionEnergy() {}

    public AlphaTorsionEnergy(TorsionList torsionList, DistanceMatrix distanceMatrix, 
		       AlphaTorsionParametersList  parametersList, double weight) {
	super(toArray(distanceMatrix, torsionList), parametersList, weight);
	this.torsionList = torsionList;
	this.distanceMatrix = distanceMatrix;
	createElementsList(torsionList);
	comment = "alphaTorsion";
    }




    public EnergyElement createElement(Object baseElement, Parameters parameters) {
	return new AlphaTorsionEnergyElement((Torsion)baseElement,
	                                    (AlphaTorsionParameters) parameters, weight);
    }

    public void handleMissingParameters(Object obj) {}


}    
