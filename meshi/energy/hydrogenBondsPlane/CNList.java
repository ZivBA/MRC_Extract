package meshi.energy.hydrogenBondsPlane;

import java.util.Iterator;

import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.MatrixRow;
import meshi.molecularElements.Atom;
import meshi.util.MeshiIterator;
import meshi.util.MeshiList;
import meshi.util.Updateable;
import meshi.util.filters.Filter;

public class CNList extends CNtwoDistancesList implements Updateable {

    //--------------------------------- data fields -----------------------------
    /*
    * The  PairsList updated every X steps.
    */
    private final int UPDATE_EVERY_X_STEPS = 50;

    protected static CNdistanceList inputNewCNList = new CNdistanceList();

    private int residuesSize;
    /*
    * holds the relevant row numbers
    */
    private MeshiList relevantRows;
    private DistanceMatrix distanceMatrix;

    /*
     * Used to avoid reupdate in the same minimization step
     */
    private int numberOfUpdates =0;
    /*
     * CNBondList/this should be updated when countUpdates >= UPDATE_EVERY_X_STEPS
     */
    private int countUpdates = 50;
    /*
     * Filter for update
     */
    private final IsGoodCNPair isGoodCNPair = new IsGoodCNPair();


    public static CNdistanceList inputNewCNList() { return inputNewCNList;}
    public final int countUpdates() {return countUpdates;}


    private CNtwoDistances [][] pairsAtResidues;
    private final IsAlive isAlive;
    //-------------------------------- constructors --------------------------------

    public CNList(DistanceMatrix distanceMatrix, int residuesSize) {
        this.distanceMatrix = distanceMatrix;
        relevantRows = new MeshiList();
        this.residuesSize = residuesSize;
        pairsAtResidues = new CNtwoDistances[residuesSize][];
        for (int i = 0; i < residuesSize; i++)
              pairsAtResidues[i] = new CNtwoDistances[i];
        isAlive =  new IsAlive(distanceMatrix );
        update_dm();

    }

    //--------------------------------------- methods ---------------------------------------

    public void update(int numberOfUpdates) {
        if (numberOfUpdates == this.numberOfUpdates+1) {
            this.numberOfUpdates++;
            update();
        }
        else if (numberOfUpdates != this.numberOfUpdates)
            throw new RuntimeException("Something weird with HbondList.update(int numberOfUpdates)\n"+
                                       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
    }


    private void update() {
        if (countUpdates == UPDATE_EVERY_X_STEPS) {
            //reset();
            countUpdates = 0;
            updateNewOnly(inputNewCNList);
            markList_dm();
            cleanList_dm();
            //dataReset();
            //update_dm();
        }
        else if (countUpdates > UPDATE_EVERY_X_STEPS)
            throw new RuntimeException("Something weird with HbondList.update()\n"+
                                       "countUpdates = "+countUpdates+" UPDATE_EVERY_X_STEPS  = "+UPDATE_EVERY_X_STEPS);
        else {
            countUpdates++;
            updateNewOnly(inputNewCNList);
            markList_dm();
        }
    }

    private void updateNewOnly(DistanceList newCNPairsList) {
        Iterator atomPairs = newCNPairsList.iterator();
        Distance dis;
        int res1, res2;
        while ((dis = (Distance) atomPairs.next()) != null) {
            res1 = dis.atom1().residueNumber();
            res2 = dis.atom2().residueNumber();
                if (pairsAtResidues[res1][res2] == null)
                        pairsAtResidues[res1][res2] = new CNtwoDistances();
                pairsAtResidues[res1][res2].setDistance(dis);
                if (pairsAtResidues[res1][res2].counter() == 2){  //only ==, if > it means pair was deleted and was formed again during UPDATE_EVERY_X_STEPS steps
                    fastAdd(pairsAtResidues[res1][res2]);
                    if (!pairsAtResidues[res1][res2].isInDM())
                        throw new RuntimeException("the alive CN pair as a not alive one");
                }
        }
    }

    public void markList_dm(){
           Iterator allPairs = iterator();
           CNtwoDistances pair;
           Distance dis1, dis2;
           boolean isDis1, isDis2;

        while((pair = (CNtwoDistances) allPairs.next()) != null){
               dis1 = pair.distance1();
               dis2 = pair.distance2();
            if (pair.isInDM()){
               isDis1 = isAlive.accept(dis1);
               isDis2 = isAlive.accept(dis2);
               if (!isDis1 || !isDis2){
                   if (!isDis1){
                        int res1 = dis1.atom1().residueNumber();
                        int res2 = dis1.atom2().residueNumber();
                        pairsAtResidues[res1][res2].deleteDistance(dis1);
                   }
                   if (!isDis2){
                        int res1 = dis2.atom1().residueNumber();
                        int res2 = dis2.atom2().residueNumber();
                        pairsAtResidues[res1][res2].deleteDistance(dis2);
                   }
               }
              }
            else{
                if (dis1 == null && dis2 != null)
                throw new RuntimeException("Wrong order of distances");
                // only both or the second distance can be == null
                if (pair.counter()!=0){
                    int res1 = dis1.atom1().residueNumber();
                    int res2 = dis1.atom2().residueNumber();
                    pairsAtResidues[res1][res2].deleteDistance(dis1);
                }
            }
        }
   }


    public void cleanList_dm(){
           Object[] newInternalArray = new Object[capacity];
           Iterator allPairs = iterator();
           CNtwoDistances pair;
           int currentIndex = 0;
           while((pair = (CNtwoDistances) allPairs.next()) != null){
               if (pair.isInDM()){
                   newInternalArray [currentIndex] = pair ;
                   currentIndex++;
               }
           }
           internalArray = newInternalArray ;
           size = currentIndex;
       }

    private void update_dm(){
        Distance dis;
        Atom headAtom;
        MatrixRow matrixRow;
        Iterator rows = distanceMatrix.rowIterator();
        Iterator rowIterator;
        CN_AtomAttribute headAttribute;
        inputNewCNList.reset();

        if (relevantRows.isEmpty()){
            // System.out.println("HBondList: relevantRows is empty");
            while((matrixRow = (MatrixRow) rows.next()) != null){
                headAtom = matrixRow.atom;
                headAttribute = (CN_AtomAttribute) headAtom.getAttribute(CN_AtomAttribute.key);
                if (headAttribute != null){//meens that it is C or N
                    relevantRows.fastAdd(new Integer(headAtom.number()));//add the row number
                    rowIterator = matrixRow.nonBondedIterator();
                    while((dis = (Distance) rowIterator.next()) != null){
                        if(isGoodCNPair.accept(dis)){
                            inputNewCNList.fastAdd(dis);
                        }
                    }//while
                }//if
            }//while there is more rows
        }//if not empty
        else {
            Iterator relevantRowsIterator = relevantRows.iterator();
            Integer current;
            int row;
            while((current = (Integer) relevantRowsIterator.next()) != null){
                row = current.intValue();
                matrixRow = distanceMatrix.rowNumber(row);
                headAtom = matrixRow.atom;
                rowIterator = matrixRow.nonBondedIterator();
                while((dis = (Distance) rowIterator.next()) != null){
                    if(isGoodCNPair.accept(dis)){
                        inputNewCNList.fastAdd(dis);
                    }
                }//while
            }//while there is more rows
        }//else
        updateNewOnly(inputNewCNList);
    }

    private void dataReset(){
        for (int i = 0; i < residuesSize; i++)
        for (int j = 0; j < i; j++)
         if (pairsAtResidues[i][j]!=null)
               pairsAtResidues[i][j]=null;
    }
 //*/
    public Iterator WithinRmaxIterator() {
        return new WithinRmaxIterator(this);// HBondIterator is an internal class.
    }


//-----------------------------------------------------------------------------------
    static class IsAlive implements Filter{
        DistanceMatrix dm;
        public IsAlive(DistanceMatrix  matrix){
            super();
            dm = matrix ;
        }

        public boolean accept(Object obj) {
            Distance distance = (Distance)obj;
            double dis = distance.distance();
            //boolean alive = distance .update(rMax2 ,rMaxPlusBuffer2 );
            if(dis == Distance.INFINITE_DISTANCE  ){  //i change it after oneStep3 - to compare it should not make any different...
             //   System.out .println("HBondList: deleted  "+obj);
                return false;
            }
            return true;
        }
    }

   //--------------------------- internal class WithinRmaxIterator ---------------------------

    private  class WithinRmaxIterator extends MeshiIterator  {
        IsWithInRmax isWithInRmax;

        public WithinRmaxIterator(MeshiList list){
            super(list);
            isWithInRmax = new IsWithInRmax();
        }
        /*
         * @return the next element that itas 2 atoms are within Rmax or Null if there is no such element
         */
        public Object next() {
            return super.next(isWithInRmax);
        }
    }


    //--------------------------- internal class IsWithInRmax ---------------------------

     class IsWithInRmax implements Filter{
        private double rMax;
        public IsWithInRmax(){
            super();
            rMax = DistanceMatrix.rMax();
        }

        public boolean accept(Object obj) {
            CNtwoDistances twoDistances = (CNtwoDistances)obj;
            if (!twoDistances.isInDM())
                return false;
            Distance dis1 = twoDistances.distance1();
            if ((rMax-dis1.distance()) >= 0){
                Distance dis2 = twoDistances.distance2();
                return ((rMax-dis2.distance()) >= 0);
            }
            else
                return false;
        }

    }

}



