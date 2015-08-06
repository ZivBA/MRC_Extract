package meshi.energy.solvate;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.hydrogenBonds.HydrogenBondDahiyatNoBBList;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/**
 * The different between this creator and the "SolvateCreatorRegularHB" is that only SC-BB or SC-SC hydrogen bonds are considered here.
 **/

public class SolvateCreatorRegularHBnoBB extends EnergyCreator  implements KeyWords {

    // The different weights relevent to this term.
    private double weightSCPolarSolvate=1.0;
    private double weightSCCarbonSolvate=1.0;
    private double weightBBPolarSolvate=1.0;
    private double weightHB=1.0;

    public SolvateCreatorRegularHBnoBB(double weightSCPolarSolvate, double weightSCCarbonSolvate, 
				       double weightBBPolarSolvate, double weightHB) {
	super(1.0);
	this.weightSCPolarSolvate = weightSCPolarSolvate;
	this.weightSCCarbonSolvate = weightSCCarbonSolvate;
	this.weightBBPolarSolvate = weightBBPolarSolvate;
	this.weightHB = weightHB;
    }
   
    public SolvateCreatorRegularHBnoBB(double weight) {
	super(weight);
	this.weightSCPolarSolvate = weight;
	this.weightSCCarbonSolvate = weight;
	this.weightBBPolarSolvate = weight;
	this.weightHB = weight;
    }

    public SolvateCreatorRegularHBnoBB() {
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

	return new SolvateEnergy(protein.atoms(), 
				 distanceMatrix, 
				 (SolvateParametersList) parametersList,
				 new HydrogenBondDahiyatNoBBList(distanceMatrix, protein.atoms(), (SolvateParametersList) parametersList),
				 weightSCPolarSolvate,
				 weightBBPolarSolvate,
				 weightSCCarbonSolvate,
				 weightHB);
    }

}
