package meshi.energy.oldSolvate;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
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

public class SolvateCreator extends EnergyCreator  implements KeyWords {

    // Default values for the HB energy term weight and the HB angular score thresholds.  
	private double simpleHBweight = 0.0;
	private double sigmoidBeginsWithH = 90.0;
	private double sigmoidEndsWithH = 100.0;
	private double sigmoidBeginsNoH = 80.0;
	private double sigmoidEndsNoH = 90.0;

    public SolvateCreator(double cooperativeSolvateWeight, 
        double simpleHBweight,
    	double sigmoidBeginsWithH, 
    	double sigmoidEndsWithH,
    	double sigmoidBeginsNoH,
    	double sigmoidEndsNoH) {
	this(cooperativeSolvateWeight,simpleHBweight);
	this.sigmoidBeginsWithH = sigmoidBeginsWithH;
	this.sigmoidEndsWithH = sigmoidEndsWithH;
	this.sigmoidBeginsNoH = sigmoidBeginsNoH;
	this.sigmoidEndsNoH = sigmoidEndsNoH;
    }

    public SolvateCreator(double cooperativeSolvateWeight, double simpleHBweight) {
	this(cooperativeSolvateWeight);
	this.simpleHBweight = simpleHBweight;
    }

    public SolvateCreator(double weight) {
	super(weight);
    }
    
    public SolvateCreator() {
	super(SOLVATE_ENERGY);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
		if (parametersList== null) {
	    	// Appanding the path to the list of parameter filenames.
	    	String[] strlist = new String[SOLVATE_PARAMETERS.length];
	    	String pathname = parametersDirectory(commands).concat("/");
	    	for (int cc=0 ; cc<SOLVATE_PARAMETERS.length ; cc++)
	        	strlist[cc] = pathname.concat(SOLVATE_PARAMETERS[cc]);
	    	parametersList = new SolvateParametersList(strlist);
	 	}
		SolvateEnergy solvateEnergy = new SolvateEnergy(protein.atoms(), 
					distanceMatrix, 
					(SolvateParametersList) parametersList,
					simpleHBweight, 
    				sigmoidBeginsWithH, 
    				sigmoidEndsWithH,
    				sigmoidBeginsNoH,
    				sigmoidEndsNoH,
					weight());
		solvateEnergy.setComment("Solvate");
		return solvateEnergy;
    }

}
