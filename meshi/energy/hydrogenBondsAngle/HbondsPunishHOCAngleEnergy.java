package meshi.energy.hydrogenBondsAngle;

import meshi.energy.hydrogenBond.HBondList;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;

/**
 * Created by IntelliJ IDEA.
 * User: amilev
 * Date: 14/03/2006
 * Time: 14:43:32
 * This class punish HydrogenBonds with HOC angles < 150 
 * The energy function is zero when the angle is >=150 or when the distance between the Hydrogen
 * and the Oxygen is bigger then 3.5 A.
 */
public class HbondsPunishHOCAngleEnergy extends AbstractPunishAngleEnergy {

    public HbondsPunishHOCAngleEnergy(){}

    public HbondsPunishHOCAngleEnergy(DistanceMatrix distanceMatrix,
                                      HBondList hBondList,
                                      double weight,
                                      DistanceList specialDis)
    {
        super(distanceMatrix,hBondList ,weight,specialDis);
    }


    /**
     * set the energyElement to point on the relavent instance
     */
    public void setEnergyElement() {
        energyElement = new HbondsPunishHOCAngleEnergyElement(distanceMatrix,weight);
    }

    /**
     * set the comment to be the relavent comment according to the instance
     */
    public void setComment() {
        comment = "HbondsPunishHOCAngleEnergy";
    }

}
