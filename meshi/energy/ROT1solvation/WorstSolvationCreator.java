package meshi.energy.ROT1solvation;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.ROT1solvation.parameters.AbstractROT1Parameters;
import meshi.energy.ROT1solvation.parameters.WorstParametersGivenCutoff;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;

public class WorstSolvationCreator extends EnergyCreator {

	private String cutoff = null;
	private boolean toCalcDerivatives = true;
	private AbstractROT1Parameters parameters = null;

	public WorstSolvationCreator(double weight, String cutoff , boolean toCalcDerivatives) {
		super(weight);
		this.toCalcDerivatives = toCalcDerivatives;
		this.cutoff = cutoff;
	}

	public WorstSolvationCreator(double weight, String cutoff) {
		this(weight, cutoff, true);
	}

	public AbstractEnergy createEnergyTerm(Protein protein,
			DistanceMatrix distanceMatrix, CommandList commands) {
		if (parameters== null) {			
	    	String pathname = parametersDirectory(commands).concat("/");
	    	String fullPathName = pathname.concat(ROT1_SOLVATE_PARAMETERS);
	    	parameters = new WorstParametersGivenCutoff(fullPathName,cutoff);
	 	}
		return new ROT1SolvationEnergy(protein.atoms(), 
				distanceMatrix,
				parameters,
				toCalcDerivatives, 
				null,
				weight());
	}

}
