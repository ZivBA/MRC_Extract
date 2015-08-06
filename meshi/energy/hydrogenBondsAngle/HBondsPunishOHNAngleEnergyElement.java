/*
 * Created on 25/01/2005
 * as part of meshi.1.5
 * 
 */
package meshi.energy.hydrogenBondsAngle;

import meshi.geometry.Angle;
import meshi.geometry.DistanceMatrix;

/**
 * @author amilev
 *
 * This class punish HydrogenBonds with angles < 150 
 * The energy function is zero when the angle is >=150 or when the distance between the Hydrogen 
 * and the Oxygen is bigger then 4 A.
 */
public class HBondsPunishOHNAngleEnergyElement extends AbstractPunishAngleEnergyElement {

    //----------------------------------- constructors --------------------------------------

    public HBondsPunishOHNAngleEnergyElement(DistanceMatrix distanceMatrix,double weight)
    {
       super(distanceMatrix, weight);
    }


    //---------------------------------- methods --------------------------------------------

    public String comment(){
        return "HBondsPunishOHNAngleEnergyElement";
    }

    public void setTheirdAtom(){
        theirdAtom = hAtom.residue().amideN();
    }

    public void setDAngleEnergyDatoms(){
                 dAngleEnergyDAngle = 2*a2*angleValue+b2;
                 dAngleEnergyDxOAtom = dAngleEnergyDAngle*angle.dangleDx1();//oAtom = atom1
                 dAngleEnergyDyOAtom = dAngleEnergyDAngle*angle.dangleDy1();
                 dAngleEnergyDzOAtom = dAngleEnergyDAngle*angle.dangleDz1();
                 dAngleEnergyDxHAtom = dAngleEnergyDAngle*angle.dangleDx2();//hatom = atom2
                 dAngleEnergyDyHAtom = dAngleEnergyDAngle*angle.dangleDy2();
                 dAngleEnergyDzHAtom = dAngleEnergyDAngle*angle.dangleDz2();
                 dAngleEnergyDxTheidAtom = dAngleEnergyDAngle*angle.dangleDx3();//natom = atom3
                 dAngleEnergyDyTheidAtom = dAngleEnergyDAngle*angle.dangleDy3();
                 dAngleEnergyDzTheidAtom = dAngleEnergyDAngle*angle.dangleDz3();
    }

    public void setAngle()
    {
              angle = new Angle(oAtom,hAtom,theirdAtom,distanceMatrix,false);//false meens don't use fakeAngle
    }
    public void updateTheirdAtom(){
       updateTheirdAtom(weight);
    }

    public void updateTheirdAtom(double weight){
        if (! theirdAtom.frozen()) {
            theirdAtom.addToFx(-1 * deDxTheidAtom * weight); // force = -derivative
            theirdAtom.addToFy(-1 * deDyTheidAtom * weight); // force = -derivative
            theirdAtom.addToFz(-1 * deDzTheidAtom * weight); // force = -derivative
        }

    }

    public String toString() {
        return "HBondsPunishOHNAngleEnergyElement (n,h,o): "+theirdAtom.residueNumber()+" "+hAtom.residueNumber()+
            " "+oAtom.residueNumber()+" dis: "+distanceValue+"\n"+
            theirdAtom+"\n"+hAtom+"\n"+oAtom+"\nangle:"+angle+"\n"
            +"disEnergy "+distanceEnergy+" angleEnergy "+angleEnergy+" energy*weight "+energy*weight;
    }

}
