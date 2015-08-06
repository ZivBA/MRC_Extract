/*
 * Created on 17/11/2004
 *
 */
package meshi.energy.hydrogenBondsPlane;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;

public class HydrogenBondsPlaneCreator extends EnergyCreator implements KeyWords {

    //------------------------- data fields ---------------------------
    public HydrogenBondsPlaneEnergy hydrogenBondsPlaneEnergy;

       /**
	 * @return Returns the hydrogenBondsEnergy.
	 */
	public final HydrogenBondsPlaneEnergy getHydrogenBondsPlaneEnergy() {	return hydrogenBondsPlaneEnergy;	}

    private String strMaxAngleOfPlainDistortion;
    private double maxAngleOfPlainDistortion;
    //------------------------- constructors ---------------------------
    
    public HydrogenBondsPlaneCreator() {
		super(HYDROGEN_BONDS_PLANE);
    }
    
	/**
	 * @param weight
	 */
	public HydrogenBondsPlaneCreator(double weight) {
		super(weight);
		key = HYDROGEN_BONDS_PLANE;
    }

    //------------------------------ methods ------------------------------
    /*
    * @param protein does not been used; must be given to implement an abstruct method of EnergyCreator.
    * @param distanceMatrix
    * @param commands gives the path to the data directory
    */
	public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
                                           CommandList commands) {
        //create parameters
		parametersList= null;
        setParametersList(commands);
        if (distanceMatrix.rMax()<5.5)
        throw new RuntimeException("Value of distanceMatrix.rMax = "+distanceMatrix.rMax()+
                " is too small for HydrogenBondsPlane." +
                "\nAdvisable value of distanceMatrix.rMax should be no smaller than 5.5");
        hydrogenBondsPlaneEnergy = new HydrogenBondsPlaneEnergy(distanceMatrix,
                                                      (HydrogenBondsPlaneParametersList) parametersList,
                                                      weight(),
                                                      maxAngleOfPlainDistortion,
                                                      new CNList(distanceMatrix, protein.residues().size()));
        distanceMatrix.energyTermsDistanceLists().fastAdd(CNList.inputNewCNList());
		return hydrogenBondsPlaneEnergy;
	}
	
    public void setParametersList(CommandList commands) {
	CommandList constrainCommands = commands.firstWordFilter(key);
	strMaxAngleOfPlainDistortion = constrainCommands.secondWord("MaxAngleOfPlaneDistortion").thirdWord();
	maxAngleOfPlainDistortion = (new Double(strMaxAngleOfPlainDistortion)).doubleValue()*Math.PI/180.0;
        if (maxAngleOfPlainDistortion >= 90)
            throw new RuntimeException("MaxAngleOfDistortion for HydrogenBondsPlane should be less than 90");

    System.out.println("MaxAngleOfDistortion = " + maxAngleOfPlainDistortion);
    }

}
