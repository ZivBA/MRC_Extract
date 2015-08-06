package meshi.geometry;
import java.util.Iterator;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.AtomPair;
import meshi.molecularElements.AtomPairList;
import meshi.molecularElements.Protein;
import meshi.util.MeshiIterator;
import meshi.util.MeshiList;
import meshi.util.Updateable;
import meshi.util.UpdateableException;
import meshi.util.filters.Filter;
/**
 **/
public class AngleList extends MeshiList implements Updateable {
    private int numberOfUpdates = 0;
    /**
     * An empty Angle list
     **/
    public AngleList() {
        this(new IsAngle());
    }
    /**
     * An empty AngleList
     **/
     protected AngleList(Filter filter) {
        super(filter);
    }
	
    public AngleList(AtomPairList bonds, DistanceMatrix distanceMatrix) {
	this();
	MeshiIterator bondIter = bonds.meshiIterator();
	AtomPair bond1, bond2;	
	while((bond1 = (AtomPair) bondIter.next()) != null) {
	    while((bond2 = (AtomPair) bondIter.nestedNextTo()) != null) 
		if (bond1.sharedAtom(bond2) != null)
		    add(new Angle(bond1,bond2,distanceMatrix));
	}
    }
    
    public void update(int numberOfUpdates) throws UpdateableException{
	if (numberOfUpdates == this.numberOfUpdates+1) {
	    int size = size();
	    for (int i = 0; i < size; i++) {
		angleAt(i).update(numberOfUpdates);
	    }
	    this.numberOfUpdates++;
	}
	else if (numberOfUpdates != this.numberOfUpdates) 
	    throw new RuntimeException("Something weird with AngleList.update(int numberOfUpdates)\n"+
				       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
    }

    public Angle angleAt(int i) { return (Angle) elementAt(i);}
    
    /**
     * Returns a sub-list containing angles that have a known name
     **/    
    public AngleList namedFilter() {
	Iterator angles = iterator();
	Angle angle;
	AngleList out = new AngleList();
	while ((angle = (Angle) angles.next()) != null) 
	    if (isNamed(angle)) out.add(angle);
	return out;
    }    
    
    public boolean isNamed(Angle angle) {
       	if (angle.getAngleName().compareTo("") != 0)
    	   return true;
    	else
    	   return false;
    }    
    
    
    public AtomList atomList() {
	AtomList list = new AtomList();
	Iterator iter = iterator();
	Angle angle;
	Atom atom1, atom2, atom3;
	while ((angle = (Angle) iter.next()) != null) {
	    atom1 = angle.atom1();
	    atom2 = angle.atom2();
	    atom3 = angle.atom3();
	    if (! list.contains(atom1)) {
		list.add(atom1);
	    }
	    if (! list.contains(atom2)) {
		list.add(atom2);
	    }
	    if (! list.contains(atom3)) {
		list.add(atom3);
	    }
	}
	return list;
    }
    public static class IsAngle implements Filter {
	public boolean accept(Object obj) {
		return (obj instanceof Angle); 
        }
    }

    public boolean sortable() {return false;}

    public static AngleList  getCaAngles(Protein protein, DistanceMatrix distanceMatrix) {
	AtomList CAs = protein.atoms().CAFilter();
	AtomPairList CaPairs = new AtomPairList();

        Atom atom1,atom2;
        MeshiIterator atomIter = CAs.meshiIterator();
        while((atom1 = (Atom) atomIter.next()) != null) {
                while((atom2 = (Atom) atomIter.nestedNextTo()) != null)
                    if (Math.abs(atom1.residueNumber() - atom2.residueNumber())<2)
			if (atom1.residueNumber() < atom2.residueNumber())
			    CaPairs.add(new AtomPair(atom1,atom2));
			else
			    CaPairs.add(new AtomPair(atom2,atom1));
        }
	
	return new AngleList(CaPairs, distanceMatrix);
    }
                                                                                                                            

    
}

