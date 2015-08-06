package meshi.energy.solvate;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.hydrogenBonds.HydrogenBondDahiyatForScmodList;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/**
 *This class creates a solvation term for the application SCMOD (concurrent sidechain modeling).
 *For this purpose the hydrogen bonds lengths are considered as long (up to 3.6Ang between donor-acceptor), 
 *and the angles of the hydrogen bonds are more relaxed (see the angle parameters in 
 *"DahiyatLowAccuracyParameters" as compared with those in "DahiyatHighAccuracyAngleParameters").
 *The derivatives and forces in this solvate implementation are never calculated. 
 *
 **/

public class SolvateCreatorSCMOD extends EnergyCreator  implements KeyWords {

	// The different weights relevent to this term.
    private double weightSCPolarSolvate=1.0;
    private double weightSCCarbonSolvate=1.0;
    private double weightBBPolarSolvate=1.0;
    private double weightHB=1.0;

    public SolvateCreatorSCMOD(double weightSCPolarSolvate, double weightSCCarbonSolvate, 
    					double weightBBPolarSolvate, double weightHB) {
		super(1.0);
		this.weightSCPolarSolvate = weightSCPolarSolvate;
		this.weightSCCarbonSolvate = weightSCCarbonSolvate;
		this.weightBBPolarSolvate = weightBBPolarSolvate;
		this.weightHB = weightHB;
    }
   
    public SolvateCreatorSCMOD(double weight) {
		super(weight);
		this.weightSCPolarSolvate = weight;
		this.weightSCCarbonSolvate = weight;
		this.weightBBPolarSolvate = weight;
		this.weightHB = weight;
    }

    public SolvateCreatorSCMOD() {
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
					new HydrogenBondDahiyatForScmodList(distanceMatrix, protein.atoms(), (SolvateParametersList) parametersList),
				    weightSCPolarSolvate,
				    weightBBPolarSolvate,
				    weightSCCarbonSolvate,
				    weightHB,
				    false /* No forces are calculated */);
    }

}
