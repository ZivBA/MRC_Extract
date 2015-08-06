package meshi.geometry;
import meshi.molecularElements.AtomPairList;
import meshi.molecularElements.Protein;
import meshi.util.MeshiIterator;
import meshi.util.MeshiList;
import meshi.util.Updateable;
import meshi.util.UpdateableException;
import meshi.util.filters.Filter;
/**
 *A list of torsion pairs, used mainly for the various two torsion energies. 
 **/
public class TorsionPairList extends MeshiList implements Updateable {
     private int numberOfUpdates = 0;
   /**
     * An empty Torsion list
     **/
    public TorsionPairList() {
        this(new IsTorsionPair());
    }
    /**
     * An empty TorsionList
     **/
     protected TorsionPairList(Filter filter) {
        super(filter);
    }


    /**
     * A TorsionPair list based on a Torsion list. Currently, only torsion pairs from the same 
     * residue are considered. Mixed torsion pairs from different residues are not treated nor 
     * created. Also, only torsions with known biological names (PHI , CHI1, etc.) are treated.
     **/
    public TorsionPairList(TorsionList torsions, DistanceMatrix distanceMatrix) {
	this();
	MeshiIterator torsionIter = torsions.meshiIterator();
	Torsion torsion1, torsion2;
	while ((torsion1 = (Torsion) torsionIter.next()) != null)
	    while ((torsion2 = (Torsion) torsionIter.nestedNextTo()) != null) {
	    	if ((torsion2.getTorsionResNum() == torsion1.getTorsionResNum()) && 
                (torsion1.getTorsionCode() > -1) && 
                (torsion2.getTorsionCode() > -1)) {
		    add(new TorsionPair(torsion1, torsion2));
		    add(new TorsionPair(torsion2, torsion1));
		}
	}
    }


    public void update(int numberOfUpdates) throws UpdateableException{
	if (numberOfUpdates == this.numberOfUpdates+1) {
	    int size = size();
	    for (int i = 0; i < size; i++) {
		torsionPairAt(i).update(numberOfUpdates);
	    }
	    this.numberOfUpdates++;
	}
	else if (numberOfUpdates != this.numberOfUpdates) 
	    throw new RuntimeException("Something weird with TorsionPairList.update(int numberOfUpdates)\n"+
				       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
    }


    public TorsionPair torsionPairAt(int i) { return (TorsionPair) elementAt(i);}
    
	
    public static class IsTorsionPair implements Filter {
	public boolean accept(Object obj) {
		return (obj instanceof TorsionPair); 
        }
    }
    
    public boolean sortable() {return false;}

    /**
     * Create a torsion pair list from a protein. Note, that many of the torsion pairs created are
     * considered not relevent, such as {PSI , CHI4}.
     **/
    public static TorsionPairList createTorsionPairList(Protein protein, DistanceMatrix distanceMatrix) {
	AtomPairList bondList = protein.bonds();
	bondList.renumber();
	AngleList angleList = new AngleList(bondList, distanceMatrix);
	TorsionList torsionList = new TorsionList(angleList, distanceMatrix);
	return new TorsionPairList(torsionList, distanceMatrix);
    }
    public static TorsionPairList createQuickAndDirtyTorsionPairList(Protein protein, DistanceMatrix distanceMatrix) {
	    AtomPairList bondList = protein.bonds();
	    bondList.renumber();
	    AngleList angleList = new AngleList(bondList, distanceMatrix);
	    TorsionList torsionList = new QuickAndDirtyTorsionList(angleList, distanceMatrix);
	    return new TorsionPairList(torsionList, distanceMatrix);
    }
}

