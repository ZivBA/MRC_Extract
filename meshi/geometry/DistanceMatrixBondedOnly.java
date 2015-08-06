package meshi.geometry;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.AtomPair;
import meshi.molecularElements.AtomPairList;
import meshi.util.UpdateableException;
/**
 
 **/
public class DistanceMatrixBondedOnly extends DistanceMatrix {
	
	
	
    private int numberOfUpdates = 0;
    private Distance[] distanceArray;

    /*-------------------------------------------- constructors -----------------------------------------------*/
    /**
     * A distance matrix with practicaly no distance cutoffs and with "hard" atoms.
     **/
   public DistanceMatrixBondedOnly(AtomList atomList, int bondedListDepth) {
       this.atomList = atomList;
       atomList.renumber();
       this.bondedListDepth = bondedListDepth;
       reset();
    }	
	
    private void reset() {
	atomArray  = atomList.toArrayOfAtoms();
	Arrays.sort(atomArray);
	bondedList = getBondedList(atomArray, bondedListDepth);
	bondedList.sort();
	distanceArray = new Distance[bondedList.size()];
	for (int i = 0; i < bondedList.size(); i++)
	    distanceArray[i] = (Distance) bondedList.elementAt(i);
	Arrays.sort(distanceArray, new DistanceComparator());
    }

    public static DistanceList getBondedList(Object[] atomArray, int depth) {
	AtomList bonded;
	Atom atom, bondedAtom;
	Iterator bondedAtoms;
	DistanceList out;
	int length = atomArray.length; 
	AtomPairList tempList = new AtomPairList();
	for (int iatom = 0; iatom < length; iatom++) {
	    atom = (Atom) atomArray[iatom];
	    bonded = getBonded(atom, depth);
	    bondedAtoms = bonded.iterator();
	    while ((bondedAtom = (Atom) bondedAtoms.next()) != null) 
		tempList.fastAdd(new AtomPair(atom, bondedAtom));
	}
	tempList.sort();
	// This instantiate the Distance objects associated with bonded atom pairs. 
	// These pairs will be ignored during the non-bonded list update.	
	AtomPair atomPair;	
	Distance distance;
	out = new DistanceList();
	Iterator atomPairs = tempList.iterator();
	while ((atomPair = (AtomPair) atomPairs.next()) != null) {
	    if (atomPair.atom1Number() > atomPair.atom2Number())
		distance = new Distance(atomPair.atom1(),atomPair.atom2());
	    else 
		distance = new Distance(atomPair.atom2(),atomPair.atom1());
	    distance.setBonded();
	    out.add(distance);
	    out.add(new DistanceMirror(distance));
	}
	return out;
    }

     public static AtomList getBonded(Atom atom, int depth) {
	AtomList out = new AtomList();
	getBonded(atom, depth, out,atom.number());
	return out;
    }
    public static void  getBonded(Atom atom, int depth, AtomList out, int rootNumber) {        
	if (depth == 0) return;
	Iterator atoms = atom.bonded().iterator();	
	Atom bondedAtom;
	while ((bondedAtom = (Atom) atoms.next()) != null) {
	    if ((rootNumber < bondedAtom.number()) &
		(! out.contains(bondedAtom)))
		out.fastAdd(bondedAtom);
	    getBonded(bondedAtom, depth-1, out, rootNumber);
	}
	
    }

    /**
     * Updates the distance matrix. 
     **/
    public void update(int numberOfUpdates) throws UpdateableException {
	if (numberOfUpdates == this.numberOfUpdates+1) {
	    this.numberOfUpdates++;
            //if (numberOfUpdates % 2 == 1)
       //    if ((numberOfUpdates < 50) ||  (numberOfUpdates % 20 == 1))
	                          update();
            
	}
	else if (numberOfUpdates != this.numberOfUpdates) 
	    throw new RuntimeException("Something weird with DistanceMatrix.update(int numberOfUpdates)\n"+
				       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
    }				       

    private void update()  throws UpdateableException{
	Iterator distances = bondedList.iterator();
	Distance distance;

	while ((distance = (Distance) distances.next()) != null) 
	    if (!distance.mirror()) distance.update();
    }
   

    /** 
     * Returns the Distance object of the parameters.
     **/
    public Distance distance(Atom atom1, Atom atom2) {
	return distance(atom1.number(), atom2.number());
    }
   public Distance distance(AtomPair atomPair) { 
	return distance(atomPair.largeNumber(), atomPair.smallNumber());
    }
    public Distance distance(int atom1Number, int atom2Number) {
	Distance distance;
	int high = distanceArray.length-1;
	int low = 0;
	int middle;
	while  (low <= high )   {
	    middle = (low+high)>>1;
	    distance =  distanceArray[middle];
	    //System.out.println("yyyyy "+atom1Number+" "+atom2Number+" | "+low+" "+high+" "+middle+" | "+distance);
	    if (distance.atom1Number > atom1Number) {high = middle - 1;}
	    else {
		if (distance.atom1Number == atom1Number) {
		    if (distance.atom2Number > atom2Number) {high = middle - 1;}
		    else {
			if (distance.atom2Number == atom2Number) return distance;
			else {low = middle + 1;}
		    }
		}
		else {low = middle + 1;}
	    }
	}
	throw new RuntimeException("distance "+atom1Number+" to "+atom2Number+" not found");
    }

    private class DistanceComparator implements Comparator {
	public int compare (Object obj1, Object obj2) {
	    Distance distance1 = (Distance) obj1;
	    Distance distance2 = (Distance) obj2;
	    if (distance1.atom1Number < distance2.atom1Number) return -1;
	    if (distance1.atom1Number > distance2.atom1Number) return 1;
	    if (distance1.atom2Number < distance2.atom2Number) return -1;
	    if (distance1.atom2Number > distance2.atom2Number) return 1;
	    return 0;
	}
    }


// 	Iterator distances = bondedList.iterator();
// 	Distance distance;
// 	while ((distance = (Distance) distances.next()) != null) 
// 	    if ((distance.atom1().number() == atom1Number) &
// 		(distance.atom2().number() == atom2Number)) return distance;
// 	throw new RuntimeException("distance "+atom1Number+" to "+atom2Number+" not found");
// 	//	return null;
//     }


    public String toString() {
	return ("DistanceMatrixBondedOnly:\n"+
		"\t number of atoms \t"+atomList.size()+
		"\t rMax \t"+rMax+
		"\t buffer\t"+buffer); 
    }


    public void testNonBondedList() {
	System.out.println("No non-bonded-list to test");
    }
}
