package meshi.energy.CAsolvate;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/**
 * The Creator class for setting up an instance of CAsolvateEnergy that gives high (but derivable) energy values
 * if in a radius 'Rad' around a certain CA atom there are more than 'MaxCANeighbors' CA atoms. If the number of
 * neighbors does not exceed 'MaxCANeighbors' the energy is 0.
 *
 * In general you need to supply the following parameters:
 * -------------------------------------------------------
 * weight
 * Rad - the radius around each CA.
 * MaxCANeighbors - The maximal number of CA neighbors in radius Rad, around each CA atoms.
 *
 * You can call different constructors, that uses different part of the default values.
 *
 **/

public class CAsolvateCreator extends EnergyCreator  implements KeyWords {

    // Default values.  
	private double Rad = 6.0; // Ang
	private double MaxCANeighbors = 5.0;

    public CAsolvateCreator(double weight, 
        double Rad,
    	double MaxCANeighbors) {
	this(weight);
	this.Rad = Rad;
	this.MaxCANeighbors = MaxCANeighbors;
    }

    public CAsolvateCreator(double weight) {
	super(weight);
    }
    
    public CAsolvateCreator() {
	super(SOLVATE_ENERGY);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
		CAsolvateEnergy caSolvateEnergy = new CAsolvateEnergy(protein.atoms(), 
					distanceMatrix, 
					null,
					Rad, 
    				MaxCANeighbors, 
					weight());
		caSolvateEnergy.setComment("CAsolvate");
		return caSolvateEnergy;
    }

}
