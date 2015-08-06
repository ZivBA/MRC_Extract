package meshi.energy;
import java.util.Iterator;

import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;

public abstract class NonBondedEnergyTerm extends AbstractEnergy{
    protected DistanceMatrix distanceMatrix;
    /**
     * Must be initialized on subclass construction.
     **/
    protected NonBondedEnergyElement energyElement=null;

    public NonBondedEnergyTerm(){super();}

    /**
     * Creates a new <code>NonBondedEnergyTerm</code> instance.
     *
     * @param updateableResources an <code>Object[]</code> value
     * @param parametersList a <code>ParametersList</code> value
     * @param weight a <code>double</code> value
     * @param distanceMatrix a <code>DistanceMatrix</code> value
     */
    public NonBondedEnergyTerm(Object[] updateableResources,
                               ParametersList  parametersList,
                               double weight,
                               DistanceMatrix distanceMatrix){
        super(updateableResources, parametersList, weight);
        this.distanceMatrix = distanceMatrix;
    }

    /**
     * Creates a new <code>NonBondedEnergyTerm</code> instance.
     *
     * @param updateableResources an <code>Object[]</code> value
     * @param weight a <code>double</code> value
     * @param distanceMatrix a <code>DistanceMatrix</code> value
     */
    public NonBondedEnergyTerm(Object[] updateableResources,
                               double weight,
                               DistanceMatrix distanceMatrix){
        super(updateableResources, weight);
        this.distanceMatrix = distanceMatrix;
    }


    public NonBondedEnergyElement getElement() {return energyElement;}
    /**
     * Evaluates energy for each distance
     * @return a sum of all energy elements
     */
    public double evaluate() {
    	if (! on) return 0.0;
    	double energy = 0;
    	DistanceList nonBondedList = distanceMatrix.nonBondedList();
    	Object[] distances = nonBondedList.internalArray();
    	int size = nonBondedList.size();
    	Distance distance;
    		
    	for (int i = 0; i < size; i++) {
    		distance = (Distance) distances[i];
    		if (!distance.frozen) {
    			energyElement.set(distance);
    			energy += energyElement.evaluate();
    		}
    	}
    	
    	return energy;
    }

    /**
     * Describe <code>evaluateAtoms</code> method here.
     */
    public void evaluateAtoms() {
    	if (on) {
    		for(Iterator nonBondedIter = distanceMatrix.nonBondedList().iterator();nonBondedIter.hasNext();){
    			Distance nonBonded = (Distance) nonBondedIter.next();
    			if (!nonBonded.frozen) {
    				energyElement.set(nonBonded);
    				energyElement.evaluateAtoms();
//    				These 3 lines are good for debugging the Excluded Volume Term
//    				if (energyElement.evaluate()>0.1)
//    					System.out.println(nonBonded.distance() + "  " + energyElement.evaluate() +  "\n"+
//    							nonBonded.atom1() + "\n" + nonBonded.atom2() + "\n");
    			}    				
    		}
    	}
    }

    /**
     * Testing of one atom in all energy elements
     * @param totalEnergy a <code>TotalEnergy</code> value
     * @param atom an criminal <code>Atom</code> value
     */
    public void test(TotalEnergy totalEnergy,Atom atom) {
	if (! on) System.out.println(""+this +" is off");
        for(Iterator nonBondedIter = distanceMatrix.nonBondedList().iterator();nonBondedIter.hasNext();){
            Distance nonBonded = (Distance)nonBondedIter.next();
	    energyElement.set(nonBonded);
	    energyElement.test(totalEnergy,atom);
    }
    }
}
