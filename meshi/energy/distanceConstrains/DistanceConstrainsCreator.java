package meshi.energy.distanceConstrains;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.Command;
import meshi.util.CommandList;
import meshi.util.KeyWords;

public class DistanceConstrainsCreator extends EnergyCreator  implements KeyWords {
    public DistanceConstrainsCreator() {
	super(DISTANCE_CONSTRAINS_ENERGY);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, CommandList commands) {
	return createEnergyTerm(protein, null, commands);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	Command command = commands.firstWordFilter(key).secondWord(INPUT_FILE);
	String distanceConstrainsFile = command.thirdWord();
	if (parametersList== null)
		parametersList = new DistanceConstrainParametersList(distanceConstrainsFile);
	return new DistanceConstrainsEnergy(protein, (DistanceConstrainParametersList) parametersList, weight());
    }
}
