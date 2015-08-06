package meshi.energy.linearRG;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/**
 *An implicit solvation energy term for all-atom models, modeling a 4.0 angs solvation shell around
 *each atom.
 **/

public class LinearRgCreator extends EnergyCreator  implements KeyWords {
        
    public LinearRgCreator(double weight) {
	super(weight);
    }

    public LinearRgCreator() {
	super(1.0);
    }


    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
		return new LinearRgEnergy(protein.atoms().filter(new AtomList.BackboneFilter()), distanceMatrix, weight());
    }

}
