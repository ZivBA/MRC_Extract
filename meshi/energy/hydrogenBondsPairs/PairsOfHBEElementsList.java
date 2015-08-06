package meshi.energy.hydrogenBondsPairs;

import java.util.Iterator;

import meshi.energy.hydrogenBond.HBondList;
import meshi.energy.hydrogenBond.HydrogenBondsEnergy;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.util.MeshiIterator;
import meshi.util.MeshiList;
import meshi.util.Updateable;
import meshi.util.filters.Filter;


/*
 * All the elements in this list are pairs of HB that are good pairs
 */
public class PairsOfHBEElementsList extends MeshiList implements Updateable {

    //--------------------------------------- data -------------------------------------------


    private static IsPairOfHydrogenBondsElements filter = new IsPairOfHydrogenBondsElements();
    /*
     * update list of the hydrogen bonds, used to update this list whenever the hBondList has been changed  
     */
    protected HBondList hBondList;
    private int numberOfUpdates;
    private boolean debug = false;

   // public final DistanceList hbeElementsList(){return hbeElementsList;}
    public final HBondList hBondList(){return hBondList;}

    private IsWithInRmax isWithinRmax = new IsWithInRmax();
    //-------------------------------------- constructors ----------------------------------------
    
    public PairsOfHBEElementsList(){
        super(filter);
    }

    public PairsOfHBEElementsList(HydrogenBondsEnergy hbe){
        super(filter);
        hBondList = hbe.hBondList();
        createPairsList(hBondList );
    }

     //-------------------------------------- methods ----------------------------------------------
    
    public void update(int numberOfUpdates) {
        if (numberOfUpdates == this.numberOfUpdates+1) {
            update();
            this.numberOfUpdates++;
        }
        else if (numberOfUpdates != this.numberOfUpdates) 
            throw new RuntimeException("Something weird with PairsOfHBEElementsList.update(int numberOfUpdates)\n"+
                                       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
    }
	    
    public void update() {
        if (hBondList.countUpdates() == 0){ //when countUpdates == 0 it meens that hBondList was reset
            if(size > 20)
                doNothing();
            updateNew(hBondList.newhBondList(),hBondList .getOldSize() );
            cleanList();
            if(debug){
                System.out.println("HBE:update:in first if:el: "+ size());
            }
        }
      else
            updateNew(hBondList.newhBondList(),hBondList .getOldSize() );
            cleanList();
    }

    private void doNothing(){} //for debug
    private void createPairsList(DistanceList hBonbds){
    	int length = hBonbds.size();
        PairOfHydrogenBondsElements setPair = new PairOfHydrogenBondsElements();
        Distance pair1 , pair2;
    	for(int i=0;i<length;i++){
            pair1 = (Distance)hBonbds.fastElementAt(i);
    		for(int j=i+1;j<length;j++)
                {    	
                    pair2 = (Distance)hBonbds.fastElementAt(j);
                    setPair.set(pair1,pair2);		
                    if (setPair.lookAtThisPair) {
                        PairOfHydrogenBondsElements pairOfHydrogenBondsElements = new PairOfHydrogenBondsElements(setPair);
                        fastAdd(pairOfHydrogenBondsElements);
                    }
                }
        }
    }
  
    public void updateNew(DistanceList newlist,int lastIndexBeforeNewElements){
        if (lastIndexBeforeNewElements == -1)
                throw new RuntimeException("problem ... this method should be called just if HBondList called to updateTwoLists !");
        if (newlist .size() !=0) {
            PairOfHydrogenBondsElements setPair = new PairOfHydrogenBondsElements();
            Distance pair1 , pair2;
            //Iterator hBondListIter = hBondList .withinRmaxIterator() ;
            for(int j=0;j<newlist.size();j++){
                for(int i=0;i<lastIndexBeforeNewElements ;i++){
                  //while((pair1 = (Distance)hBondListIter .next()) !=null){
                    pair1 = (Distance) hBondList .fastElementAt(i);
                    pair2 = (Distance)(newlist.fastElementAt(j));
                    setPair.set(pair1, pair2);
                    if(setPair.lookAtThisPair) {
                        PairOfHydrogenBondsElements pairOfHydrogenBondsElements = new PairOfHydrogenBondsElements(setPair);
                        fastAdd(pairOfHydrogenBondsElements);
                    }
                }//while
            }//for
            createPairsList(newlist);
        }//if
    }

     /**
        * save only the live elements inthis list
        */
       public void cleanList(){
           Object[] newInternalArray = new Object[capacity ];
           Iterator allPairs = iterator() ;
           PairOfHydrogenBondsElements pair;
           int currentIndex =0;
           while((pair = (PairOfHydrogenBondsElements) allPairs.next() ) != null){
                if (isWithinRmax .accept(pair )){
                    newInternalArray [currentIndex] = pair ;
                   currentIndex++;
                }
              // else
//                System.out.println("PairsOfHBEElementsList: delete "+pair );
           }
           internalArray = newInternalArray ;
           size = currentIndex;
       }

    public Iterator withinRmaxPairsIterator() {
        return new WithinRmaxPairsIterator(this);
    } 

    //---------------------------------------------------------------------------------------
    
    //----------------------------- parivate class IsPairOfHydrogenBondsElements ------------

    //---------------------------------------------------------------------------------------
        private static class IsPairOfHydrogenBondsElements implements Filter {
        public boolean accept(Object obj) {
            boolean ans = ((obj instanceof PairOfHydrogenBondsElements) &&
                           ((PairOfHydrogenBondsElements)obj).lookAtThisPair);
            return ans;
        }
    }


    //--------------------------- internal class IsWithInRmax ---------------------------

    static class IsWithInRmax implements Filter{
        private double rMax;
        public IsWithInRmax(){
            super();
            rMax = DistanceMatrix.rMax() ;
        }

        public boolean accept(Object obj) {
            PairOfHydrogenBondsElements pair = (PairOfHydrogenBondsElements)obj;
            if ((rMax >= pair.HOelement1.distance()) &
                    (rMax >= pair.HOelement2.distance()))
            return true;
            else
            return false;
        }
    }
    
    //---------------------------------------------------------------------------------------
    
    //----------------------------- parivate class WithinRmaxPairsIterator ----------- ------------

    //---------------------------------------------------------------------------------------
    private  class WithinRmaxPairsIterator extends MeshiIterator  {

        private IsWithInRmax isWithInRmax;
        public WithinRmaxPairsIterator(MeshiList list) {
            super(list);
            isWithInRmax = new IsWithInRmax();
        }
        
        /*
         * @return the next element that in each of its 2 pairs the 2 atoms are within Rmax or Null if there is no such element
         */
        public Object next() {
            return next(isWithInRmax);
        }
    }
}
