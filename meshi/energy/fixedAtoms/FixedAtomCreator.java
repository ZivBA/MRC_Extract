package meshi.energy.fixedAtoms;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/**
 **/

public class FixedAtomCreator extends EnergyCreator  implements KeyWords {

    AtomList atoms;

    public FixedAtomCreator(AtomList al , double weight) {
	super(weight);
	atoms = al;
    }

    public FixedAtomCreator(AtomList al) {
	super(1.0);
	atoms = al;
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
		FixedAtomEnergy fixedAtomEnergy = new FixedAtomEnergy(atoms, distanceMatrix, weight());
		return fixedAtomEnergy;
    }

}
