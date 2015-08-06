package meshi.energy.solvate;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.hydrogenBonds.HydrogenBondDahiyatList;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/**
 * The Creator class for setting up an instance of SolvateEnergy as described in Kalisman & Keasar (2006).
 *
 * In general you need to supply the following parameters:
 * -------------------------------------------------------
 * weight - The CooperativeEnergyTerm part weight. Note, that weight has no affect on the regular HB part, that
 * has its own weight.
 * simpleHBweight - The weight of the HB part.
 * sigmoidBeginsWithH,sigmoidEndsWithH - The transition angles (in degrees) for the HB sigmoid when the hydrogen 
 * in the HB is defined in MESHI. Above sigmoidEndsWithH the sigmoid is given a value of 1.0 . Bellow          
 * sigmoidBeginsWithH the sigmoid is given a value of 0.0 . In between it raises smoothly by cubic spline.
 * sigmoidBeginsNoH,sigmoidEndsNoH - The same as above, only for HB sigmoids where the hydrogen in the HB is present.
 *
 *
 * You can call different constructors, that uses different part of the default values.
 *
 **/

public class SolvateCreatorAlt2 extends EnergyCreator  implements KeyWords {

	// The different weights relevent to this term.
    private double weightSCPolarSolvate=1.0;
    private double weightSCCarbonSolvate=1.0;
    private double weightBBPolarSolvate=1.0;
    private double weightHB=1.0;

    public SolvateCreatorAlt2(double weightSCPolarSolvate, double weightSCCarbonSolvate, 
    					double weightBBPolarSolvate, double weightHB) {
		super(1.0);
		this.weightSCPolarSolvate = weightSCPolarSolvate;
		this.weightSCCarbonSolvate = weightSCCarbonSolvate;
		this.weightBBPolarSolvate = weightBBPolarSolvate;
		this.weightHB = weightHB;
    }
   
    public SolvateCreatorAlt2(double weight) {
		super(weight);
		this.weightSCPolarSolvate = weight;
		this.weightSCCarbonSolvate = weight;
		this.weightBBPolarSolvate = weight;
		this.weightHB = weight;
    }

    public SolvateCreatorAlt2() {
		super(SOLVATE_ENERGY);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
		if (parametersList== null) {
	    	// Appanding the path to the list of parameter filenames.
	    	String[] strlist = new String[SOLVATE_PARAMETERS_ALT2.length];
	    	String pathname = parametersDirectory(commands).concat("/");
	    	for (int cc=0 ; cc<SOLVATE_PARAMETERS_ALT2.length ; cc++)
	        	strlist[cc] = pathname.concat(SOLVATE_PARAMETERS_ALT2[cc]);
	    	parametersList = new SolvateParametersList(strlist);
	 	}

		return new SolvateEnergy(protein.atoms(), 
					distanceMatrix, 
					(SolvateParametersList) parametersList,
					new HydrogenBondDahiyatList(distanceMatrix, protein.atoms(), (SolvateParametersList) parametersList),
				    weightSCPolarSolvate,
				    weightBBPolarSolvate,
				    weightSCCarbonSolvate,
				    weightHB);
    }

}
