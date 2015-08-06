package meshi.energy.localRG;
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

public class LocalRgCreator extends EnergyCreator  implements KeyWords {
	
	private int begins=-1,ends=-1;
        
    public LocalRgCreator(double weight , int begins , int ends) {
	super(weight);
	this.begins = begins;
	this.ends = ends;
    }

    public LocalRgCreator() {
	super(1.0);
    }


    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
		return new LocalRgEnergy(protein.atoms().filter(new AtomList.BackboneFilter()), distanceMatrix, begins, ends, weight());
    }

}
