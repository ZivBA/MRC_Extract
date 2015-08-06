package meshi.molecularElements;
import java.util.Iterator;

import meshi.util.SortableMeshiList;
import meshi.util.filters.Filter;

public class AtomPairList extends SortableMeshiList {
    /**
     * An empty AtomPair list
     **/
    public AtomPairList() {
	this(new IsAtomPair());
    }
    /**
     * An empty AtomPair list
     **/
    protected AtomPairList(Filter filter) {
	super(filter);
    }
    public AtomPairList(ResidueList residueList) {
	this();
	Iterator residues, bonds;
	Residue residue;
	AtomPair bond;
	Residue prevResidue = null;
	residues = residueList.iterator();
	while ((residue = (Residue)residues.next()) != null) {
	    bonds = residue.bonds().iterator();
	    while ((bond = (AtomPair) bonds.next()) != null)
		add(bond);
	    if (prevResidue != null) {
		if (prevResidue.nextAtom() != null) {  // Current and previous residues must be connected. 
		    if ((prevResidue.nextAtom() != residue.tail()) |
			(prevResidue.head() != residue.prevAtom())) {
			if ((prevResidue.head() != null) &
			    (residue.tail() != null)) {
			    add(prevResidue.head().bond(residue.tail()));
			    prevResidue.setNextAtom(residue.tail());
			    residue.setPrevAtom(prevResidue.head());
			}
		    }
		    else add(new AtomPair(residue.tail(),prevResidue.head()));
		}
		else {
		    if ((prevResidue.head() != null) &
			(residue.tail() != null)) {
			add(prevResidue.head().bond(residue.tail()));
			prevResidue.setNextAtom(residue.tail());
			residue.setPrevAtom(prevResidue.head());
		    }
		}
	    }
	    prevResidue = residue;
	}
	renumber();
    }
	    
    public AtomPair atomPairAt(int i) { 
	return (AtomPair) elementAt(i);
    }
    
          // fastAdd is equal to add(Object element) in MeshiList , but a checking of type is excluded
    protected boolean fastAdd(AtomPair element) {
	if (size < capacity) {
	    internalArray[size] = element;
	    size++;
	    modCount++;
	    return true;
	}
	else {
	    capacity *= 2;
	    Object[] newArray = new Object[capacity];
	    for (int i = 0; i < size; i++)
		newArray[i] = internalArray[i];
	    internalArray = newArray;
	    return fastAdd(element) ;
	}
    }
    
    public AtomList atomList() {
	AtomList list = new AtomList();
	Iterator iter = iterator();
	AtomPair atomPair;
	Atom atom1, atom2;
	while ((atomPair = (AtomPair) iter.next()) != null) {
	    atom1 = atomPair.atom1();
	    atom2 = atomPair.atom2();
	    if (! list.contains(atom1)) {
		list.add(atom1);
	    }
	    if (! list.contains(atom2)) {
		list.add(atom2);
	    }
	}
	return list;
    }
    public static class IsAtomPair implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof AtomPair); 
	}
    }
    public boolean sortable() {return true;}

    public void renumber() {
	Iterator iter = iterator();
	AtomPair atomPair; 
	while ((atomPair = (AtomPair) iter.next()) != null) 
	    atomPair.renumber();
    }

    public void printLaconic(String prompt) {
	Iterator iter = iterator();
	AtomPair atomPair; 
	while ((atomPair = (AtomPair) iter.next()) != null) 
	    System.out.println(atomPair.laconic(prompt));
    }	
}
