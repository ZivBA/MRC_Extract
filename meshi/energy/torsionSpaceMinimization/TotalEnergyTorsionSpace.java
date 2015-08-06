package meshi.energy.torsionSpaceMinimization;

import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;

public class TotalEnergyTorsionSpace extends TotalEnergy {
	
	private TorsionMinimizationTree tmt;
	
    public TotalEnergyTorsionSpace(Protein protein, 
            DistanceMatrix distanceMatrix,
            EnergyCreator[] energyCreators,
            CommandList commands) {
    	super(protein, distanceMatrix, energyCreators, commands);
    	if (frozenAtomsExist())
    		throw new RuntimeException("Frozen atoms exists in the protein. Torsion space cannot work with frozen atoms.");
    	tmt = new TorsionMinimizationTree(protein, distanceMatrix);
    	setCoordinates();		
    }

	
    public double evaluate() {
    	resetAtomForces();
    	tmt.buildProtein();
    	super.evaluate();
    	tmt.calcDerivatives();
    	return totalEnergy; 
    }
     
	protected void setCoordinates(AtomList atomList) {
		coordinates = null;
	}

	protected void setCoordinates() {
		coordinates = tmt.getCoors();
	}	
	
	public void resetAtomForces() {
		for (int c=0 ; c<atomList.size(); c++)
			atomList.atomAt(c).coordinates().resetForces();
	}
	
}
