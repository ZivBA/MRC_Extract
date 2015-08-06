package meshi.molecularElements.hydrogenBonds;
import java.util.Iterator;
import java.util.Vector;

import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.Updateable;
import meshi.util.UpdateableException;
import meshi.util.filters.KolDichfin;


/** 
 * A super class with the most general treatment of a list of hydrogen bonds in an atom list. 
 * A specific implementation should put content to these abstract methods according to the specific
 * implementaion of the hydrogen bonds. The constructor gets a distance matrix object, and an 
 * AtomList object of all the atoms where bonds are searched. The class is centered around the
 * "bondList" vector that stores all the hydrogen bond objects relevent to the atom list. Whenever 
 * a new hydrogen bond appears in the non-bonded-list, it is added to this vector. Hydrogen bonds that are
 * no longer in the non-bonded list are occasionaly removed from this vector. The parameter 
 * 'refreshVectorEvery' in the constructor determine how often this 'clean-up' is done. The bondList
 *  vector is non-redundent.
 * 
 *
 * IPORTANT IMPORTANT IMPORTANT
 * NOTE: The distance matrix must be in an updated state before any of this class's
 * methods are called. If it is not up to date, the results would be meaningless!!!
 * This class does not update the distance matrix on its own at any time!!
 *
 **/
public abstract class AbstractHydrogenBondList implements Updateable {

	// Fields:
	// -------
	protected int maxAtomNum;
	protected int atomListSize;
    protected DistanceMatrix dm=null; 
    protected AtomList atomList=null; 
	private DistanceList newToNonBonded; // A pointer to this list is added to the DistanceMatrix's list of lists (see line 74). This way, only new distances are recognized. 
    protected Vector<AbstractHydrogenBond> bondList;
    protected int[] lut; 	// The look-up table (lut) converts the atom internal number (field of Atom) to its index in the atom list given to the constructor.      
    protected AbstractHydrogenBond[][] lutHB; // lutHB[i] is an array of pointers to all the hydrogen bond objects where atom number i is either their donor or acceptor. It is null for non-polar atoms, or atoms not participating in a hydrogen bond.     
    private int numberOfUpdates = 0; // for the Updateable interface
    private int refreshVectorEvery=200;


    public AbstractHydrogenBondList() {}

    public AbstractHydrogenBondList(DistanceMatrix dm, AtomList atomList,int refreshVectorEvery) {
    	this.dm = dm;
    	this.atomList = atomList;
    	this.refreshVectorEvery = refreshVectorEvery;
    	bondList = new Vector<AbstractHydrogenBond>();

		// Building the look-up tables
		maxAtomNum=-1;
		atomListSize = atomList.size();
		for (int c=0; c<atomListSize ; c++) {
		    if (atomList.atomAt(c).number() > maxAtomNum)
				maxAtomNum = atomList.atomAt(c).number();
		}
		lut = new int[maxAtomNum+1];
		lutHB = new AbstractHydrogenBond[maxAtomNum+1][];
		for (int c=0; c<maxAtomNum ; c++) {
		    lut[c] = -1;
		}
		for (int c=0; c<atomListSize ; c++) {
		    lut[atomList.atomAt(c).number()] = c;
		}		
		
		// This part is needed so that we can get the new distances that were added
		// to the non-bonded list in the last update.
		newToNonBonded = new DistanceList(new KolDichfin());
		dm.energyTermsDistanceLists().fastAdd(newToNonBonded);
		
		buildSpecificStructures();
		
    }

    
    public void update(int updateNumber) throws UpdateableException {
        if (updateNumber == numberOfUpdates+1) {
        	if ((updateNumber % refreshVectorEvery == 0) && (updateNumber>0)) 
        		update(true,updateNumber);
        	else
        		update(false,updateNumber);
            this.numberOfUpdates++;
        }
        else if (updateNumber != this.numberOfUpdates)
            throw new RuntimeException("Something weird with the update in AbstractHydrogenBondList.java\n"+
                    "The update number is "+updateNumber+" while the numberOfUpdates in the class is "+numberOfUpdates);
    }    

    
    /** 
     * The 'toRefreshVector' parameter determines if broken H-bonds are removed from the vector.
     **/
    protected void update(boolean toRefreshVector, int updateNumber) throws UpdateableException {
    	addNewBondsToList();
        for(AbstractHydrogenBond hb : bondList)
            hb.updateHBvalueAndDerivatives();
    	if (toRefreshVector)
    		removeBrokenHB();
    }
    

	/**
	* This methods add new h-bonds to the bondList vector. Care is taken to keep the vector non-redundent.
	**/ 
	private void addNewBondsToList() {
		Iterator iter;
		if (numberOfUpdates == 0)   // The first update to the H-bond list
			iter = dm.nonBondedList().iterator();
		else
			iter = newToNonBonded.iterator();
		Distance dis; 
		AbstractHydrogenBond hb;
		while((dis = (Distance) iter.next()) != null) {
	    	// Possible HB Ahoy!! (atom1 and atom2 are capable of forming a hydrogen bond) 
    		if (!dis.atom1().isCarbon && !dis.atom2().isCarbon && !dis.atom1().isHydrogen && !dis.atom2().isHydrogen) {
   	    		// We first need to check that this pair of polars is not already in the vector.
   	    		if (findBondByPolars(dis.atom1(),dis.atom2())==null) {
   	    			hb = createHBfromPolars(dis.atom1(),dis.atom2());
   	    			if (hb!=null) {
   	    				bondList.add(hb);
   	    				addToLUTvec(lutHB,dis.atom1(),hb);
   	    				addToLUTvec(lutHB,dis.atom2(),hb);
   	    			}
   	    		}
   	    	}
   	    }
   	    newToNonBonded.reset();
	}


	/** 
	 * This method remove from the vector Hydrogen Bond objects, whose distance between polars is
	 * so large, they no longer appears in the non-bonded list of the distance matrix.
	 **/ 
	private void removeBrokenHB(){
		Vector<AbstractHydrogenBond> tmpVector = new Vector<AbstractHydrogenBond>();
		AbstractHydrogenBond[][] tmpLutHB = new AbstractHydrogenBond[lutHB.length][]; 
		AbstractHydrogenBond hb;
		DistanceList dislist = dm.nonBondedList();
		Iterator iter = dislist.iterator();
		Distance dis; 
		while((dis = (Distance) iter.next()) != null) {
	    	// Possible HB Ahoy!! (atom1 and atom2 are capable of forming a hydrogen bond) 
    		if (!dis.atom1().isCarbon && !dis.atom2().isCarbon && !dis.atom1().isHydrogen && !dis.atom2().isHydrogen) {
    	    		hb = findBondByPolars(dis.atom1(),dis.atom2());
    	    		if (hb!=null) {
    	    			tmpVector.add(hb);
    	    			addToLUTvec(tmpLutHB,dis.atom1(),hb);
    	    			addToLUTvec(tmpLutHB,dis.atom2(),hb);
    	    		}
    	    }
    	}
    	bondList.clear();
    	bondList = tmpVector;
    	lutHB = tmpLutHB;
	}

	/** 
	 * This method returns a pointer to the hydrogen bond that exists between two POLAR atoms (hydrogens are not consider polar)
	 * If no such object exists then null is returned.
	 **/ 
	public AbstractHydrogenBond findBondByPolars(Atom atom1 , Atom atom2) {
   	   	if (lutHB[atom1.number()]!=null)
   	    	for (int c=0 ; c<lutHB[atom1.number()].length ; c++)
   	    		if ((lutHB[atom1.number()][c].getFirstPolar().number() == atom2.number()) ||
   	    			(lutHB[atom1.number()][c].getSecondPolar().number() == atom2.number()))
   	    			return lutHB[atom1.number()][c];
   	    return null;
	}


	/** 
	 *Auxilarry method: Updating a vector, 'vec' similar to the lutHB field, with the new HB 
	 *The reason we do not access lutHB directly, is that this method is also called from the 
	 *'removeBrokenHB' method, and there a different vector need updating.
	 **/ 
	private void addToLUTvec(AbstractHydrogenBond[][] vec, Atom atom, AbstractHydrogenBond hb) {
		if (vec[atom.number()] == null) {
			vec[atom.number()] = new AbstractHydrogenBond[1];
			vec[atom.number()][0] = hb;
		}
		else {
			AbstractHydrogenBond[] tmp = new AbstractHydrogenBond[vec[atom.number()].length+1];
			for (int c=0 ; c<vec[atom.number()].length ; c++)
				tmp[c] = vec[atom.number()][c];
			tmp[vec[atom.number()].length] = hb;
			vec[atom.number()] = tmp;	
		}
	}

	public void print() {
		for (int c=0 ; c<bondList.size() ; c++) {
			System.out.println(bondList.get(c));
		}
	}


    /** 
      * This method should build any implementation-specific data structures needed.
     **/
    protected abstract void buildSpecificStructures();
    	
    /** 
      * This method should build an implementaion specific HB from two polar atoms given as parameters.
      * If the creation is not possible then null is returned.
     **/
    protected abstract AbstractHydrogenBond createHBfromPolars(Atom atom1,Atom atom2);

	
}


  

