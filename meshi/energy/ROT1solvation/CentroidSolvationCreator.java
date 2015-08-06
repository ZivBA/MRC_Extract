package meshi.energy.ROT1solvation;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.ROT1solvation.parameters.AbstractROT1Parameters;
import meshi.energy.ROT1solvation.parameters.CentroidParametersGivenCutoff;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;

public class CentroidSolvationCreator extends EnergyCreator {

	private String cutoff = null;
	private boolean toCalcDerivatives = true;
	private AbstractROT1Parameters parameters = null;
	private boolean[] matchingRes = null;

	public CentroidSolvationCreator(double weight, String cutoff , boolean toCalcDerivatives){
	this(weight, cutoff ,toCalcDerivatives, null);
}
	
	public CentroidSolvationCreator(double weight, String cutoff , boolean toCalcDerivatives, boolean[] matchingRes) {
		super(weight);
		this.toCalcDerivatives = toCalcDerivatives;
		this.cutoff = cutoff;
		this.matchingRes = matchingRes;
	}

	public CentroidSolvationCreator(double weight, String cutoff) {
		this(weight, cutoff, true,null);
	}

	public void setMatchingResidues(boolean[] newMatcing){
		matchingRes = newMatcing;
	}

	public AbstractEnergy createEnergyTerm(Protein protein,
			DistanceMatrix distanceMatrix, CommandList commands) {
		if (parameters== null) {			
	    	String pathname = parametersDirectory(commands).concat("/");
	    	String fullPathName = pathname.concat(ROT1_SOLVATE_PARAMETERS);
	    	parameters = new CentroidParametersGivenCutoff(fullPathName,cutoff);
	 	}
		return new ROT1SolvationEnergy(protein.atoms(), 
				distanceMatrix,
				parameters,
				toCalcDerivatives, 
				matchingRes,
				weight());
	}

}