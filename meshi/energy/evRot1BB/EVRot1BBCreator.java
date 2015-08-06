package meshi.energy.evRot1BB;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;


public class EVRot1BBCreator extends EnergyCreator  implements KeyWords {

	private double[][] pp = null;
	private double maximalZ=1e10;

/*    public EVRot1BBCreator(double weight, double maximalZ) {
        super(weight);
        this.maximalZ = maximalZ;
    }


    public EVRot1BBCreator(double maximalZ) {
        super(SOLVATE_ENERGY);
        this.maximalZ = maximalZ;
    }
*/
    public EVRot1BBCreator(double weight,double[][] pp , double maximalZ) {
	super(weight);
	this.pp = pp;
	this.maximalZ = maximalZ;
    }


    public EVRot1BBCreator(double[][] pp , double maximalZ) {
	super(SOLVATE_ENERGY);
	this.pp = pp;
	this.maximalZ = maximalZ;
    }
    
    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
/*
if (pp==null) {
	DunbrackLib lib = new DunbrackLib(commands, 1.0 , 2);
	pp = RotamericTools.putIntoRot1(protein, distanceMatrix, lib);
}
*/

                if (parametersList== null) {
	    	// Appanding the path to the list of parameter filenames.
	    	String pathname = parametersDirectory(commands).concat("/");
	    	parametersList = new EVRot1BBParameterList(pathname.concat(EXCLUDED_VOL_PARAMETERS));
	 	}
		EVRot1BBEnergy evEnergy = new EVRot1BBEnergy(protein.atoms(), 
					distanceMatrix, 
					(EVRot1BBParameterList) parametersList,
					pp,
					maximalZ,
					weight());
		evEnergy.setComment("Rot1 EV BB");
		return evEnergy;
    }

}
