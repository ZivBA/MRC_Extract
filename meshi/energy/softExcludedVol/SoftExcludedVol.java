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
package meshi.energy.softExcludedVol;
import meshi.energy.NonBondedEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.util.Classes;
public class SoftExcludedVol extends NonBondedEnergyTerm implements Classes{
    public SoftExcludedVol(){super();}

    public SoftExcludedVol(DistanceMatrix distanceMatrix, 
                           SoftExcludedVolParametersList  parametersList, 
                           int type,
                           double weight){
	super(toArray(distanceMatrix), parametersList, weight,distanceMatrix);
	comment = "SoftExcludedVol";
	energyElement = new SoftExcludedVolEnergyElement(parametersList, distanceMatrix, type, weight);
    }
}

	
 
