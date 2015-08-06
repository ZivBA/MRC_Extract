package meshi.energy.hydrogenBondsPlane;

import java.util.Iterator;

import meshi.energy.NonBondedEnergyTerm;
import meshi.energy.TotalEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.parameters.AtomTypes;

public class HydrogenBondsPlaneEnergy extends NonBondedEnergyTerm implements AtomTypes {

    //--------------------------------- data fields ----------------------------------
    protected CNList CNList; //element list
    private HydrogenBondsPlaneEnergyElement energyElement;
    private HydrogenBondsPlaneParametersList parametersList;

    public CNList CNList() {
        return CNList;
    }
    /**
     * Returns the parametersList.
     */
    public final HydrogenBondsPlaneParametersList getParametersList() {
        return parametersList;
    }

    //------------------------------------ constractor ---------------------------------------------

    public HydrogenBondsPlaneEnergy (DistanceMatrix distanceMatrix,
                                     HydrogenBondsPlaneParametersList parametersList,
                                     double weight,
                                     double maxAngleOfDistortion,
                                     CNList CNList){
        super(toArray(distanceMatrix,CNList), //create the array of apdateableResources
              parametersList,
              weight,
              distanceMatrix);
        comment = "HB-PlaneEnergy";
        this.distanceMatrix = distanceMatrix;
        this.CNList = CNList;
        energyElement = new HydrogenBondsPlaneEnergyElement(distanceMatrix, weight, maxAngleOfDistortion);
        this.parametersList = parametersList;
        Iterator atoms = distanceMatrix.atomList().iterator();
        Atom atom;
        while ((atom = (Atom) atoms.next()) != null) {
            if (IsCN.isC(atom)) atom.addAttribute(new CN_AtomAttribute(true,false));
            if (IsCN.isN(atom)) atom.addAttribute(new CN_AtomAttribute(false,true));
        }
    }


    //-------------------------------------- methods ---------------------------------------

    public double evaluate() {
        if (! on) return 0.0;
        double energy = 0;
        Iterator cnIter = CNList.WithinRmaxIterator();
        CNtwoDistances pairDis;
        while ((pairDis  = (CNtwoDistances) cnIter.next()) != null) {
               energyElement.set(pairDis);
               energy += energyElement.evaluate();
           }
        if ((! (energy < 0) ) & (! ( energy == 0)) & (!(energy > 0)))
            System.out.println("weird energy "+energy);

         return energy;
    }

    public void evaluateAtoms() {
        if (! on) return;
        Iterator cnIter = CNList.WithinRmaxIterator();
        CNtwoDistances pairDis;
        while ((pairDis  = (CNtwoDistances) cnIter.next()) != null) {
               energyElement.set(pairDis);
               energyElement.evaluate();
           }
    }

    public void test(TotalEnergy totalEnergy, Atom atom) {
        if (! on) System.out.println(""+this +" is off");
        Iterator cnIter = CNList.WithinRmaxIterator();
        CNtwoDistances pairDis;
        while ((pairDis  = (CNtwoDistances) cnIter.next()) != null) {
                energyElement.set(pairDis);
                energyElement.setAtoms();
                energyElement.test(totalEnergy,atom);
            }    
    }
}
