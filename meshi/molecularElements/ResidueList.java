package meshi.molecularElements;
import java.util.Iterator;

import meshi.parameters.Residues;
import meshi.sequences.SequenceAlignmentCell;
import meshi.util.SortableMeshiList;
import meshi.util.filters.Filter;

/**
 * A list of residues.
 * Note that some of the residues may be dummy. 
 * Specifically, the first residue may be dummy in order to be compatible with the 
 * biologists convention that the first residue is 1 (they grew up on FORTRAN).  
 **/
public class ResidueList extends SortableMeshiList{
    /**
     * The residue number of the first residue (may but also may not be 1).
     **/
    protected int firstNonDummyResidueIndex;
    protected Residue firstNonDummyResidue;

    // --------------------------------- constructors ------------------------------
    /**
     * Empty list.
     **/    public ResidueList() {
	super(new IsResidue());
	firstNonDummyResidueIndex = -1;
	firstNonDummyResidue = null;
    }

    /**
     * Constructs a ResidueList from a sequence (a string of one letter codes).
     * The residueCreator is responsible for the interpretation of the letters to actual residues and thus
     * determines the molecular model. 
     **/
    public ResidueList(String sequence, ResidueCreator creator) {
	this();
	add(new DummyResidue(0)); //traditionally sequences start at 1
	char[] seq = sequence.toCharArray();
	if (seq.length < 1) throw new RuntimeException("Protein sequence length needs to be at least one");
	if (seq.length < 2)
	    add(creator.create(String.valueOf(seq[0]),1,Residues.SINGLE, 0.0, 0.0, 0.0));
	else {
	    add(creator.create(String.valueOf(seq[0]),1,Residues.NTER, 0.0, 0.0, 0.0));
	    for (int i = 1; i < seq.length-1; i++)
		add(creator.create(String.valueOf(seq[i]),i+1,Residues.NORMAL, 4.0*i, 0.0, 0.0));
	    add(creator.create(String.valueOf(seq[seq.length-1]),seq.length,Residues.CTER, 
			       (seq.length-1)*10.0, 0.0, 0.0));
	}
	sort();
    }

   /**
     * Constructs a ResidueList from a list of atoms.
     * The residueCreator allows the use of information that is not  stored in the atoms themselves. 
     **/
    public ResidueList(AtomList atomList, ResidueCreator creator) {
	this();
	Atom atom;
	Iterator atoms = atomList.iterator();
	AtomList newAtomList = null;
	boolean first = true;
	add(new DummyResidue(0)); //traditionally sequences start at 1
	if (! atoms.hasNext()) throw new RuntimeException(" No Atoms in AtomList "+atomList.comment());
	while ((atom = (Atom) atoms.next ()) != null) {
	    while (size()< atom.residueNumber()) {
		if (newAtomList == null) add(new DummyResidue(size()));
		else {
		    if (first) {
			add(creator.create(newAtomList,size(),Residues.NTER));
			first = false;
		    }
		    else add(creator.create(newAtomList,size(),Residues.NORMAL));
		    newAtomList = null;
		}
	    }
	    if (size() == atom.residueNumber()) {
		if (newAtomList == null) newAtomList = new AtomList();
		newAtomList.add(atom);
	    }
	}
	add(creator.create(newAtomList,size(),Residues.CTER));
	sort();
     }

		
    //------------------------------------------- add -----------------------------------
    public boolean add(Object element) {
	Residue residue = (Residue) element;
	if (! residue.dummy()) {
	    if ((firstNonDummyResidue == null) ||
		(residue.number < firstNonDummyResidueIndex)){
		firstNonDummyResidue = residue;
		firstNonDummyResidueIndex = residue.number;
	    }
	}
	return(super.add(element));
    }

    //------------------------------------------- methods -----------------------------------
    /**
     * Fetch a residue by its position in the list.
     **/
    public Residue residueAt(int index) {
	return (Residue) elementAt(index);
    }
    
    /**
     * Fetch a residue by its residue number.
     **/
     public Residue residue(int residueNumber) {
	Residue key = new Residue(residueNumber);
	int index = binarySearch(key);
	if (index < 0) return null;
	return residueAt(index);
    }
	
    /**
     * Extracts the residues accepted by the filter.
     **/
    public ResidueList filter(Filter residueFilter) {
	Iterator iter = iterator();
	Residue residue;
	ResidueList newList = new ResidueList();
	while((residue = (Residue) iter.next()) != null) {
	    if (residueFilter.accept(residue))
		newList.add(residue);
	}
	return newList;
    }

    static class IsResidue implements Filter {
	public boolean accept(Object obj) {return (obj instanceof Residue);}
    } 


    public int firstNonDummyResidueIndex() {return firstNonDummyResidueIndex;}
    public Residue firstNonDummyResidue() {return firstNonDummyResidue;}

    public  int firstNonDummyResidueNumber() {
	int ires;
	int size = size();
	for (ires = 0; ires < size; ires++)
	    if (! (elementAt(ires) instanceof DummyResidue)) return ires;
	return -1;	
    }
    
    public int numberOfNonDummyResidues() {
	int out = 0;
	for (Iterator residues = iterator(); residues.hasNext();) {
	    Residue residue = (Residue) residues.next();
	    if (!residue.dummy()) out++;
	}
	return out;
    }

	    

    public void initiateAllResidues() {
        Iterator residueIter = iterator();
        Residue residue;
        while ((residue=(Residue)residueIter.next())!=null){
            residue.initiateAtoms();
        }
    }

    public AtomList atoms() {return new AtomList(this);}

    public String toString() {
	String out = "";
	Iterator residues = iterator(); 
	Residue residue = (Residue) residues.next(); // remove the first dummy residue
	while (residues.hasNext()) {
	    residue = (Residue) residues.next();
	    if (residue.dummy()) out+=SequenceAlignmentCell.WILDCARD_CHAR;
	    else out+=residue.nameOneLetter();
	}
	return out;
    }
}
