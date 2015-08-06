package meshi.energy.hydrogenBond;

import java.util.Iterator;

import meshi.energy.NonBondedEnergyTerm;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.parameters.AtomTypes;

/**
 * @author amilev
 *
 * This class claculate the HydrogenBonds Energy over all the unbounded Hydrogen-Oxygen Pairs in the protein. 
 * To do so we keep an updatable-pairs-list with all the unbounded Hydrogen-Oxygen Pairs.
 */
public class HydrogenBondsEnergy extends NonBondedEnergyTerm implements AtomTypes {

    //--------------------------------- data fields ----------------------------------

    //*
    // * Accept a Distance between atoms if one of them is backbond Hydrogen and the other is backbond Oxygen
    // */
    //protected DistanceMatrix distanceMatrix;
    //protected Distance elementsList;
    protected HBondList hBondList;
    protected HydrogenBondsEnergyElement energyElement;
    private HydrogenBondsParametersList parametersList;
   // private DistanceList specialDisatnces = null;
    /**
     * @return Returns the distanceMatrix.
     */
    public final DistanceMatrix getDistanceMatrix() {
        return distanceMatrix;
    }

    public HBondList hBondList() {
        return hBondList;
    }
    /**
     * @return Returns the parametersList.
     */
    public final HydrogenBondsParametersList getParametersList() {
        return parametersList;
    }

    //------------------------------------ constractor ---------------------------------------------

    public HydrogenBondsEnergy(){super();}


    public HydrogenBondsEnergy (DistanceMatrix distanceMatrix,
                                HydrogenBondsParametersList parametersList,
                                double weight,
                                HBondList hBondList){
        super(toArray(distanceMatrix,hBondList), //create the array of apdateableResources
              parametersList,
              weight,
              distanceMatrix);
        comment = "HydrogenBondsEnergy";
        this.distanceMatrix = distanceMatrix;
        this.hBondList = hBondList;
        energyElement = new HydrogenBondsEnergyElement(weight);
        this.parametersList = parametersList;
    }

   public HydrogenBondsEnergy (DistanceMatrix distanceMatrix,
                                HydrogenBondsParametersList parametersList,
                                double weight,
                                HBondList hBondList,
                                DistanceList specialDistances){
       this(distanceMatrix,parametersList,weight ,hBondList );
    //   this.specialDisatnces = specialDistances ;

   }


    //-------------------------------------- methods ---------------------------------------

    public double evaluate() {
        if (! on) return 0.0;
        double energy = 0, energyCurrent;
        Iterator hoIter = hBondList.withinRmaxIterator();
        Distance pair;

      /*  if(specialDisatnces != null){
            Iterator specDisIter = specialDisatnces.iterator() ;
            while((pair  = (Distance) specDisIter .next()) != null){
                 if(DistanceMatrix.rMax() > pair .distance() ){
                HB_DistanceAttribute distanceAttribute = new HB_DistanceAttribute(true);
                distanceAttribute.setParameters((HydrogenBondsParameters) parametersList.parameters(pair));
                distanceAttribute .set(pair.atom1() ,pair .atom2() );
                pair.addAttribute(distanceAttribute);
                    energyElement .set(pair );
                    energyCurrent = energyElement.evaluate();
                    energy += energyCurrent;
                    energyElement.freeElement();
                 if(! (energy >0 | energy <= 0) )
            System.out.print("eee");
                 }

            }
        }*/

        while ((pair  = (Distance) hoIter.next()) != null) {
            energyElement.set(pair);
            energyCurrent = energyElement.evaluate();
            energy += energyCurrent;
            energyElement.freeElement();
             if(! (energy >0 | energy <= 0) )
            System.out.print("eee");
        }

         if(! (energy >0 | energy <= 0) )
            System.out.print("eee");
         return energy;

    }

    public void evaluateAtoms() {
       if (! on) return;
        Iterator hoIter = hBondList.withinRmaxIterator();
        Distance pair;
        while ((pair  = (Distance) hoIter.next()) != null) {
            energyElement.set(pair);
            energyElement.evaluateAtoms();
            energyElement.freeElement();
        }
    }
}
