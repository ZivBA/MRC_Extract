/**
 *Excluded Volume potential.
 *
 *A potential that creates a strong repulsion between atoms when they are getting
 *near their VDW radius. There is no attraction part like in the VDW.
 *The functional form of the term is:  
 *
 *dis =                EV
 *
 *[0,sigma]  C*(dis-sigma)^4
 *[sigma,Inf]          0
 *
 *ALPHA is set in ExcludedVolParameters.java. Currently it is 0.2 Ang.
 *ALPHA is the transition zone (in Ang) where the energy change in the forth power 0.0 to 1.0.
 *
 **/
package meshi.energy.excludedVol;
import java.util.Iterator;

import meshi.energy.NonBondedEnergyTerm;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.util.Classes;
import meshi.util.filters.Filter;
public class ExcludedVol extends NonBondedEnergyTerm implements Classes{
    protected double Rfac;
    protected Filter filter; 

    public ExcludedVol(){super();}

    public ExcludedVol(DistanceMatrix distanceMatrix, 
                       ExcludedVolParametersList  parametersList, 
                       double Rfac , double weight, Filter filter){
	super(toArray(distanceMatrix), parametersList, weight,distanceMatrix);
	comment = "ExcludedVol";
	energyElement = new ExcludedVolEnergyElement(parametersList, distanceMatrix, Rfac, 
                                                     weight);
	this.filter = filter;
    }

    public double evaluate() {
	if (! on) return 0.0;
	double e, energy = 0;
	DistanceList nonBondedList;
	if (filter == null) 
		nonBondedList = distanceMatrix.nonBondedList();
	else {
	    nonBondedList = new DistanceList();
	    distanceMatrix.nonBondedList().filter(filter,nonBondedList);
	}
	Object[] distances = nonBondedList.internalArray();
	int size = nonBondedList.size();
	Distance distance;
	for (int i = 0; i < size; i++) {
	    distance = (Distance) distances[i];
	    energyElement.set(distance);
	    energy += energyElement.evaluate();
        }
	return energy;
    }
    
    public void evaluateAtoms() {
	if (on) {
		Iterator distanceListIter = distanceMatrix.nonBondedList().iterator();
		Distance distance;
		while ((distance  = (Distance) distanceListIter.next()) != null) {
	    	energyElement.set(distance);
			energyElement.evaluateAtoms();
	    }
	}
    }
}

	
 
