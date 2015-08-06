package meshi.energy.solvateNoHB;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;


public class SolvateNoHBCreator extends EnergyCreator  implements KeyWords {



    public SolvateNoHBCreator(double weight) {
	super(weight);
    }
    
    public SolvateNoHBCreator() {
	super(SOLVATE_ENERGY);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
		if (parametersList== null) {
	    	// Appanding the path to the list of parameter filenames.
	    	String[] strlist = new String[SOLVATE_NOHB_PARAMETERS.length];
	    	String pathname = parametersDirectory(commands).concat("/");
	    	for (int cc=0 ; cc<SOLVATE_NOHB_PARAMETERS.length ; cc++)
	        	strlist[cc] = pathname.concat(SOLVATE_NOHB_PARAMETERS[cc]);
	    	parametersList = new SolvateNoHBParametersList(strlist);
	 	}
		SolvateNoHBEnergy solvateEnergy = new SolvateNoHBEnergy(protein.atoms(), 
					distanceMatrix, 
					(SolvateNoHBParametersList) parametersList,
					weight());
		solvateEnergy.setComment("NoHB Solvate");
		return solvateEnergy;
    }

}
