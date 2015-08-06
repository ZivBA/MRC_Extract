package meshi.energy.rot1Pairwise;

import meshi.energy.NonBondedEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.util.Classes;

public class Rot1Pairwise extends NonBondedEnergyTerm implements Classes{
    public Rot1Pairwise(){super();}

    public Rot1Pairwise(DistanceMatrix distanceMatrix, 
            Rot1PairwiseParametersList  parametersList, 
            int type,
            double weight){
	super(toArray(distanceMatrix), parametersList, weight,distanceMatrix);
	comment = "Rot1Pairwise";
	this.distanceMatrix = distanceMatrix;
	energyElement = new Rot1PairwiseEnergyElement(parametersList, type, weight);
    }
}
