/*
 * Created on 26/01/2005
 * as part of meshi.1.5
 * 
 */
package meshi.energy.hydrogenBondsAngle;

import meshi.energy.hydrogenBond.HBondList;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;

/**
 * @author amilev
 *
 * This class punish HydrogenBonds with angles < 150 
 * The energy function is zero when the angle is >=150 or when the distance between the Hydrogen 
 * and the Oxygen is bigger then 3.5 A.
 */
public class HBondsPunishOHNAngleEnergy extends AbstractPunishAngleEnergy {
	
    public HBondsPunishOHNAngleEnergy(){}

    public HBondsPunishOHNAngleEnergy(DistanceMatrix distanceMatrix,
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
        energyElement = new HBondsPunishOHNAngleEnergyElement(distanceMatrix,weight);
    }

    /**
     * set the comment to be the relavent comment according to the instance
     */
    public void setComment() {
        comment = "HBondsPunishOHNAngleEnergy";
    }

}
