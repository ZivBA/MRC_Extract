/*
 * Created on 26/12/2004
 *
 */
package meshi.energy.hydrogenBondsPairs;

import java.util.Iterator;

import meshi.energy.NonBondedEnergyTerm;
import meshi.energy.TotalEnergy;
import meshi.energy.hydrogenBond.HB_DistanceAttribute;
import meshi.energy.hydrogenBond.SynchronizedUpdateException;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.parameters.AtomTypes;

/**
 * @author amilev
 */
public class HydrogenBondsPairsEnergy extends NonBondedEnergyTerm implements AtomTypes {

    //--------------------------------- data ------------------------------------------

    public PairsOfHBEElementsList getPairsOfHBEElementsList() {
        return pairsOfHBEElementsList;
    }

    //protected HydrogenBondsPairsEnergyElement energyElement;
    protected PairsOfHBEElementsList pairsOfHBEElementsList;


    private static final double DIFULT_PUNISHMENT = -10;        //TODO check this parameter
    private static final double DIFULT_HPUNISHMENT = -100;    //TODO check this parameter

    //--------------------------------- constructors ----------------------------------
    
    /**
     * 
     */
    public HydrogenBondsPairsEnergy() {
        super();	
    }

    public HydrogenBondsPairsEnergy(DistanceMatrix distanceMatrix,
                                    HelixParametersList helixParametersList,
                                    BetaParametersList betaParametersList,
                                    PairsOfHBEElementsList pairsOfHBEElementsList,
                                    double weight) {
        this(distanceMatrix,helixParametersList,betaParametersList, pairsOfHBEElementsList, weight,DIFULT_PUNISHMENT,DIFULT_HPUNISHMENT);
    }

    public HydrogenBondsPairsEnergy(DistanceMatrix distanceMatrix,
                                    HelixParametersList helixParametersList,
                                    BetaParametersList betaParametersList,
                                    PairsOfHBEElementsList pairsOfHBEElementsList,
                                    double weight,
                                    int [] specialDisArray,
                                    boolean antiParallel) {
        this(distanceMatrix,helixParametersList,betaParametersList, pairsOfHBEElementsList, weight,DIFULT_PUNISHMENT,DIFULT_HPUNISHMENT,specialDisArray,antiParallel);
    }

    public HydrogenBondsPairsEnergy(DistanceMatrix distanceMatrix,
                                    HelixParametersList helixParametersList,
                                    BetaParametersList betaParametersList,
                                    PairsOfHBEElementsList pairsOfHBEElementsList,
                                    double weight,
                                    double punish,
                                    double hpunish ){

        super(toArray(distanceMatrix,pairsOfHBEElementsList.hBondList(),pairsOfHBEElementsList),
              weight,
              distanceMatrix);
        comment = "HB-PairsEnergy";
        this.pairsOfHBEElementsList = pairsOfHBEElementsList;
        energyElement = new HydrogenBondsPairsEnergyElement(helixParametersList,betaParametersList, weight,punish,hpunish);
    }

    public HydrogenBondsPairsEnergy(DistanceMatrix distanceMatrix,
                                    HelixParametersList helixParametersList,
                                    BetaParametersList betaParametersList,
                                    PairsOfHBEElementsList pairsOfHBEElementsList, 
                                    double weight,
                                    double punish,
                                    double hpunish,
                                    int[] specialDisArray,
                                    boolean antiParallel){
        
        super(toArray(distanceMatrix,pairsOfHBEElementsList.hBondList(),pairsOfHBEElementsList),
              weight,
              distanceMatrix);
        comment = "HB-PairsEnergy";
        this.pairsOfHBEElementsList = pairsOfHBEElementsList;
        energyElement = new HydrogenBondsPairsEnergyElement(helixParametersList,betaParametersList, weight,punish,hpunish,specialDisArray,antiParallel);
    }
	

    //---------------------------------------- methods ------------------------------------------
    
    /* (non-Javadoc)
     * @see meshi.energy.AbstractEnergy#evaluate()
     */
    public double evaluate() {
        if (! on) return 0.0;
        double energy = 0;
        Iterator pairsIter = pairsOfHBEElementsList.withinRmaxPairsIterator();
        PairOfHydrogenBondsElements pair;

        while ((pair  = (PairOfHydrogenBondsElements) pairsIter.next()) != null) {
            //if (pair.atomPair1.withinRmax() && pair.atomPair2.withinRmax())
            if (pair.isHelixPair() | pair.isBetaPair()) //TODO change when fix the all problem
            {
                energyElement.set(pair);
                energy += energyElement.evaluate();
                ((HydrogenBondsPairsEnergyElement )energyElement) .freeElenet();
            }
        }
        return energy;
    }

    /* (non-Javadoc)
     * @see meshi.energy.AbstractEnergy#evaluateAtoms()
     */
    public void evaluateAtoms() {
        if (! on) return;
        Iterator pairsIter = pairsOfHBEElementsList.withinRmaxPairsIterator();
        PairOfHydrogenBondsElements pair;
        while ((pair  = (PairOfHydrogenBondsElements) pairsIter.next()) != null) {
            {
                 if (pair.isHelixPair() | pair.isBetaPair()){ //TODO change when fix the all problem
                       energyElement.set(pair);
                       energyElement.evaluateAtoms();
                       ((HydrogenBondsPairsEnergyElement )energyElement) .freeElenet();
                 }
            }
        }
    }

    /* (non-Javadoc)
     * @see meshi.energy.AbstractEnergy#test(meshi.energy.TotalEnergy, boolean, boolean)
     */
    public void test(TotalEnergy totalEnergy,Atom atom){
        //System.out.println("Start Test "+comment);
         System.out.println("hBondList: ");
        for(int i =0;i<pairsOfHBEElementsList .hBondList .size();i++){
                    Distance currentDis = (Distance)pairsOfHBEElementsList.hBondList .elementAt(i);
                    if (currentDis.isDistanceOf(atom))
                        System.out.println(currentDis);
        }

        System.out.println("pairsOfHBEElementsList: ");
        PairOfHydrogenBondsElements pair;
        for(Iterator pairsIter = pairsOfHBEElementsList.withinRmaxPairsIterator();pairsIter.hasNext();){
            pair = (PairOfHydrogenBondsElements)pairsIter.next();
            if(pair.HOelement1.isDistanceOf(atom) | pair.HOelement2 .isDistanceOf(atom) )
                        System.out.println(pair);
        }

        if (! on) System.out.println(" "+this +" is off");
        Iterator pairsIter = pairsOfHBEElementsList.withinRmaxPairsIterator();

       while((pair = (PairOfHydrogenBondsElements)pairsIter.next()) != null){
             energyElement.set(pair);
            if(((HydrogenBondsPairsEnergyElement )energyElement).distanceAttributes1 != ((HydrogenBondsPairsEnergyElement )energyElement).HOelement1.getAttribute(HB_DistanceAttribute.key))
                System.out.println("this is wird");
            if(energyElement .atoms().whereIs(atom) >= 0){
                System.out.println(comment+" element-pair: ");
                energyElement.atoms().print();

                //System.out.println(energyElement.oAtom1());
                if(pair.lookAtThisPair){
                    //System.out.println(comment+": lookAtThisPair");
                   //System.out.println("atoms: "+ pair.atoms());
                   //System.out.println("size "+pair.atoms().size());
                   //System.out.println("pairs "+pair.atoms());
                }
                else{
                    System.out.println("dont look at this pair" +comment);
                    //pair.atoms().print();
                }
                //System.out.println("start Test pair !!!!!!!!!!!!!!!!!!");
                try{
                    energyElement.test(totalEnergy,atom);
                }
                catch(SynchronizedUpdateException e)   {
                    System.out.println(e);

                }
                //System.out.println("End Test pair !!!!!!!!!!!!!!!!!!!!!");

            }
           ((HydrogenBondsPairsEnergyElement )energyElement).freeElenet();
        }
        System.out.println("End Test "+comment);
    }
    /*
      OLD TEST 
    public void test(TotalEnergy totalEnergy, boolean testAll, boolean verbose) {
        double e;
        pairsOfHBEElementsList.print();
    	if (! on) System.out.println(""+this +" is off");
    	double highestEnergy = -100000;
    	double lowestEnergy = 100000;
    	Iterator pairsIter = pairsOfHBEElementsList.iterator();
    	PairOfHydrogenBondsElements pair;
    	while ((pair  =  (PairOfHydrogenBondsElements)pairsIter.next()) != null) {
            if (!(DistanceMatrix.rMax()-pair.HOelement1.distance() >=0) || !(DistanceMatrix.rMax()-pair.HOelement2.distance() >=0)) continue;
            System.out.println("PairOfHydrogenBondsElements:"+pair);
            energyElement.set(pair);
            e = energyElement.evaluate();
            if (e > highestEnergy) {
                System.out.println("********* Highest energy so far "+e+" ************");
                highestEnergy = e;
                energyElement.test(totalEnergy, verbose);	
            }
            else if(testAll) energyElement.test(totalEnergy, verbose);
            if (e < lowestEnergy) {
                System.out.println("********* Lowest energy so far "+e+" ************");
                lowestEnergy = e;
            }    		
    	}
    }
    */
}
