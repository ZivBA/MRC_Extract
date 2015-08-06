package meshi.energy.hydrogenBond;

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

/**
 * @author amilev
 *   Note: if the 2 atoms of the HB are frozen - we don't add them to the list.
 **/
public class HBondList extends HBdistanceList implements Updateable {

    //--------------------------------- data fields -----------------------------

    
    /*
     * The  PairsList updated every X steps.
     */
    private final int UPDATE_EVERY_X_STEPS = 50;
    ///*
    // * List of the candidate to hydrogen bonding given by the DistanceMatrix object
    // */
    //protected DistanceList nonBondedList;
    /*
     * List of all the new HB elements that were added to the  hBondList in the last update call.
     */
    protected DistanceList newhBondList;
    protected static HBdistanceList inputNewHBList = new HBdistanceList();
    /*
    * List of all the HB elements
    */
    //protected DistanceList  hBondList;

    /*
     * holds the relevant row numbers
     */
    private MeshiList relevantRows;
    /*
     * get the parameters (epsilon, sigma) of each element
     */
    private HydrogenBondsParametersList parametersList;
    private DistanceMatrix distanceMatrix;

    /*
     * Used to avoid reupdate in the same minimization step
     */
    private int numberOfUpdates =0;
    /*
     *  hBondList/this should be updated when countUpdates >= UPDATE_EVERY_X_STEPS
     */
    private int countUpdates = 50;
    /*
     * Filter for update
     */
    private GoodResiduesForHB goodResiduesForHB;
    private final IsAlive isAlive;

    /**
     * to find easily the last element that was in this list before the new update
     */
    private int oldSize = -1;

   public final int getOldSize() { return oldSize; }
    public final DistanceList newhBondList() { return newhBondList;}
    public static HBdistanceList inputNewHBList() { return inputNewHBList;}
    public final int countUpdates() {return countUpdates;}

    //-------------------------------- constructors --------------------------------

    /**
	 * @param distanceMatrix
	 */
	public HBondList(DistanceMatrix distanceMatrix,
                     HydrogenBondsParametersList parametersList) {
        super(new GoodResiduesForHB());
        goodResiduesForHB = new GoodResiduesForHB();
        this.parametersList = parametersList;
        this.distanceMatrix = distanceMatrix;
        //hBondList = new DistanceList();
        newhBondList = new DistanceList();
        relevantRows = new MeshiList();
        //nonBondedList = distanceMatrix.nonBondedList();
        //update(nonBondedList);
        isAlive =  new IsAlive(distanceMatrix );
        update_dm();
    }

    public HBondList(DistanceMatrix distanceMatrix,
                     HydrogenBondsParametersList parametersList,
                     DistanceList specialDistances) {
        super(new GoodResiduesForHB(specialDistances));
        goodResiduesForHB = new GoodResiduesForHB(specialDistances);
        this.parametersList = parametersList;
        this.distanceMatrix = distanceMatrix;
        newhBondList = new DistanceList();
        relevantRows = new MeshiList();
        isAlive =  new IsAlive(distanceMatrix );
        update_dm();
    }

    //--------------------------------------- methods ---------------------------------------

	/* (non-Javadoc)
	 * @see meshi.util.Updateable#update(int)
	 */
	public void update(int numberOfUpdates) {
		if (numberOfUpdates == this.numberOfUpdates+1) {
		    this.numberOfUpdates++;
            //System.out.println("HBondList: in update(int numberOfUpdates):");
		    update();
		}
		else if (numberOfUpdates != this.numberOfUpdates)
		    throw new RuntimeException("Something weird with HbondList.update(int numberOfUpdates)\n"+
                                       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
	}


    private void update() {
        oldSize = -1;
        newhBondList.reset();
        if (countUpdates == UPDATE_EVERY_X_STEPS) {//TODO can delete this if - i leave it for now for possible future needs
            int prevrousSize = size() ;
            updateTwoLists(inputNewHBList);
            int secondSize = size() ;
            countUpdates = 0;
            cleanList_dm();
            oldSize = size - inputNewHBList .size() ;
            //System.out.println("HBondList: prevouse size "+prevrousSize +" secondSize  "+secondSize +" current size "+size);
            if(size > secondSize | secondSize < prevrousSize | size < oldSize )
                throw new RuntimeException("HBondsList: problem at update - the list after adding elements is shorter or after cleaning is longer");
        }
        else if (countUpdates > UPDATE_EVERY_X_STEPS)
            throw new RuntimeException("Something weird with HbondList.update()\n"+
                                       "countUpdates = "+countUpdates+" UPDATE_EVERY_X_STEPS  = "+UPDATE_EVERY_X_STEPS);
        else {
            //System.out.println("HBondList: in update: distanceMatrix.newHydrogenBondsList():");
            //distanceMatrix.newHydrogenBondsList().print();
            countUpdates++;
           updateTwoLists(inputNewHBList);
           cleanList_dm();
           oldSize = size - inputNewHBList .size() ;
           if(size < oldSize)
                     throw new RuntimeException("HBondList: size < old size");
            }

    }

	/**
	 * @param atomPairList of new distances that where added after  HbondList was updated in the last time
	 */
    private void updateTwoLists(HBdistanceList atomPairList) {
        oldSize = -1;
        Iterator atomPairs = atomPairList.iterator();
	    Distance pair;
        oldSize = size;
        while ((pair = (Distance) atomPairs.next()) != null ){
            if (!isAlive.accept(pair))
                  System.out.println("HBondList: there is a problem");
            //System.out.println("HBondList: in updateTwoLists: pair: "+pair);
            //if (goodResiduesForHB.accept(pair)) {     this term is checked on HBdistanceList
            //System.out.println("!!!!!! HBondList: in updateTwoLists: Goodpair !!!!! : "+pair);
                HB_DistanceAttribute distanceAttribute = new HB_DistanceAttribute(true);
                distanceAttribute.setParameters((HydrogenBondsParameters) parametersList.parameters(pair));
                distanceAttribute .set(pair.atom1() ,pair .atom2() );
                pair.addAttribute(distanceAttribute);
                fastAdd(pair);
                newhBondList.fastAdd(pair);
             //}
        }
    }

    /*
     * the update is done by first go over the heads atoms of the row in distance matrix
     * and just then (only if the head atom is Hydrogen or Oxygen) go over the rows themself.
     * (more afficient ?)
     */
    private void update_dm(){
        // System.out.println("HBondList: in update_dm: start");
        Distance pair;
        Atom headAtom;
        MatrixRow matrixRow;
        Iterator rows = distanceMatrix.rowIterator();
        Iterator rowIterator;
        HB_AtomAttribute headAttribute;
        if (relevantRows.isEmpty()){
            // System.out.println("HBondList: relevantRows is empty");
            while((matrixRow = (MatrixRow) rows.next()) != null){
                headAtom = matrixRow.atom;
                headAttribute = (HB_AtomAttribute) headAtom.getAttribute(HB_AtomAttribute.key);
                if (headAttribute != null){//meens that it is H or O !
                    // System.out.println("HBondList: headAttribute != null");
                    relevantRows.fastAdd(new Integer(headAtom.number()));//add the row number
                    rowIterator = matrixRow.nonBondedIterator();
                    while((pair = (Distance) rowIterator.next()) != null){
                        if(goodResiduesForHB.accept(pair) ){
                            if(!isAlive .accept(pair ))
                                throw new RuntimeException("need the second check !");
                            HB_DistanceAttribute distanceAttribute = new HB_DistanceAttribute(true);
                            distanceAttribute.setParameters
                                ((HydrogenBondsParameters) parametersList.parameters(pair));
                            distanceAttribute .set(pair.atom1() ,pair .atom2() );
                            pair.addAttribute(distanceAttribute);
                            //hBondList.fastAdd(pair);
                            fastAdd(pair);
                        }
                    }//while
                }//if
            }//while there is more rows
        }//if not empty
        else {
            // System.out.println("HBondList: relevantRows is NOT empty");
            Iterator relevantRowsIterator = relevantRows.iterator();
            Integer current;
            int row;
            while((current = (Integer) relevantRowsIterator.next()) != null){
                row = current.intValue();
                matrixRow = distanceMatrix.rowNumber(row);
                headAtom = matrixRow.atom;
                rowIterator = matrixRow.nonBondedIterator();
                while((pair = (Distance) rowIterator.next()) != null){
                    if(goodResiduesForHB.accept(pair)){
                         if(!isAlive .accept(pair ))
                                throw new RuntimeException("need the secod check !");
                        HB_DistanceAttribute distanceAttribute = new HB_DistanceAttribute(true);
                        distanceAttribute.setParameters
                            ((HydrogenBondsParameters) parametersList.parameters(pair));
                        pair.addAttribute(distanceAttribute);
                        // hBondList.fastAdd(pair);
                        fastAdd(pair);
                    }
                }//while
            }//while there is more rows
        }//else
        // System.out.println("HBondList: in update_dm: end");
        // this.print();
    }


    public Iterator newWithinRmaxIterator() {
        return new WithinRmaxIterator(newhBondList);// WithinRmaxIterator is an internal class.
    }

     public boolean selectionAdd(Object object){
      return selectionAdd(object,new IsWithInRmax());
     }
    /**
        * save only the live elements inthis list
        */
       public void cleanList_dm(){
           Object[] newInternalArray = new Object[capacity ];
           Iterator allPairs = iterator() ;
           Distance pair;
           int currentIndex =0;
           while((pair = (Distance) allPairs.next() ) != null){
               if (isAlive .accept(pair )){
                   newInternalArray [currentIndex] = pair ;
                   currentIndex++;
               }
           }
           internalArray = newInternalArray ;
           size = currentIndex;
       }

       public Iterator withinRmaxIterator() {
           return new WithinRmaxIterator(this);// WithinRmaxIterator is an internal class.
       }


    //--------------------------- internal class IsWithInRmax ---------------------------

    static class IsWithInRmax implements Filter{
       private final double rMax;
        public IsWithInRmax(){
            super();
            rMax =   DistanceMatrix.rMax();
        }

        public boolean accept(Object obj) {
            Distance distance = (Distance)obj;
            //double dis = rMax - distance.distance();
            if (rMax >= distance.distance() )
                return true;
            else
                return false;
        }
    } //--------------------------- internal class IsAlive ---------------------------

    static class IsAlive implements Filter{
        DistanceMatrix dm;
        boolean firstTimeWornning = true;
        public IsAlive(DistanceMatrix  matrix){
            super();
            dm = matrix ;
        }

        public boolean accept(Object obj) {
            Distance distance = (Distance)obj;
            boolean frozen = distance.atom1() .frozen() & distance.atom2() .frozen() ;
            double dis = distance.distance();
            //boolean alive = distance .update(rMax2 ,rMaxPlusBuffer2 );
            if (frozen){
                if(firstTimeWornning )   //TODO check this - maybe frozen atoms should be looked at ?
                {
                    System.out.println("*** NOTE: frozen atoms does not get into HBondList !!! *****");
                    firstTimeWornning = false;
                }

                return false;
            }
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
}




