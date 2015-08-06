/**
 *Excluded Volume potential, somewhat following (in its repulsive part) that of Summa and Levitt (2007)
 *
 **/
package meshi.energy.excludedVolumeImprovedDistance;
import meshi.energy.NonBondedEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.util.Classes;
public class EVenergy extends NonBondedEnergyTerm implements Classes{
    public EVenergy(){super();}

    public EVenergy(DistanceMatrix distanceMatrix, 
                           EVenergyParametersList  parametersList, 
                           int type,
                           double weight){
	super(toArray(distanceMatrix), parametersList, weight,distanceMatrix);
	comment = "ExcludedVol";
	energyElement = new EVenergyElement(parametersList, distanceMatrix, type, weight);
    }
}

	
 
