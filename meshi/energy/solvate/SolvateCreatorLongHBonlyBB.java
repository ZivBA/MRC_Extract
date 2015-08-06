package meshi.energy.solvate;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.hydrogenBonds.HydrogenBondDahiyatJustBBList;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/**
 * This method is the same as 'SolvateCreatorLongHB` except that only backbone atoms are considered.   
 **/

public class SolvateCreatorLongHBonlyBB extends EnergyCreator  implements KeyWords {

	// The different weights relevent to this term.
    private double weightSCPolarSolvate=1.0;
    private double weightSCCarbonSolvate=1.0;
    private double weightBBPolarSolvate=1.0;
    private double weightHB=1.0;

    public SolvateCreatorLongHBonlyBB(double weightSCPolarSolvate, double weightSCCarbonSolvate, 
    					double weightBBPolarSolvate, double weightHB) {
		super(1.0);
		this.weightSCPolarSolvate = weightSCPolarSolvate;
		this.weightSCCarbonSolvate = weightSCCarbonSolvate;
		this.weightBBPolarSolvate = weightBBPolarSolvate;
		this.weightHB = weightHB;
    }
   
    public SolvateCreatorLongHBonlyBB(double weight) {
		super(weight);
		this.weightSCPolarSolvate = weight;
		this.weightSCCarbonSolvate = weight;
		this.weightBBPolarSolvate = weight;
		this.weightHB = weight;
    }

    public SolvateCreatorLongHBonlyBB() {
		super(SOLVATE_ENERGY);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
		if (parametersList== null) {
	    	// Appanding the path to the list of parameter filenames.
	    	String[] strlist = new String[SOLVATE_LONG_HB_PARAMETERS.length];
	    	String pathname = parametersDirectory(commands).concat("/");
	    	for (int cc=0 ; cc<SOLVATE_LONG_HB_PARAMETERS.length ; cc++)
	        	strlist[cc] = pathname.concat(SOLVATE_LONG_HB_PARAMETERS[cc]);
	    	parametersList = new SolvateParametersList(strlist);
	 	}
 	
		return new SolvateEnergy(protein.atoms(), 
					distanceMatrix, 
					(SolvateParametersList) parametersList,
					new HydrogenBondDahiyatJustBBList(distanceMatrix, protein.atoms(), (SolvateParametersList) parametersList),
				    weightSCPolarSolvate,
				    weightBBPolarSolvate,
				    weightSCCarbonSolvate,
				    weightHB);
    }

}
