package meshi.sequences;
import java.util.Iterator;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.Rms;
import meshi.util.filters.Filter;
import meshi.util.filters.KolDichfin;

public class AtomAlignment extends Alignment{
    public AtomAlignment() {
	super(new IsAlignmentColumn());
    }
    public AtomAlignment(AtomList atomList0,AtomList atomList1) {
	super(new IsAlignmentColumn());
	if (atomList0.size() != atomList1.size()) {
//		atomList0.print();
//		System.out.println("-----------------------------\n-----------------------------");
//		atomList1.print();
//		System.out.println(atomList0.size() + " " + atomList1.size());
		throw new RuntimeException("List lengths must be equal");
	}
	Iterator atoms1 = atomList1.iterator();
	for (Iterator atoms0 = atomList0.iterator(); atoms0.hasNext();)
	    add(new AtomAlignmentColumn((Atom) atoms0.next(), atomList0.comment(),
					(Atom) atoms1.next(), atomList1.comment()));
    }

    public AtomAlignment(ResidueAlignment residueAlignment) {
	    this(residueAlignment, new KolDichfin());
    }

    public AtomAlignment(ResidueAlignment residueAlignment, Filter filter) {
	super(new IsAlignmentColumn());
	for(Iterator columns = residueAlignment.iterator(); columns.hasNext();) {
	    ResidueAlignmentColumn column = (ResidueAlignmentColumn) columns.next();
	    if (!column.hasGap()) {
		    Atom ca0 = column.residue0().ca();
		    Atom ca1 = column.residue1().ca();
		    if (filter.accept(ca0) & filter.accept(ca1)) {
		    	add(new AtomAlignmentColumn(ca0, ca1));
		    }
	    }
	}
    }
		

    private static class IsAlignmentColumn implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof AtomAlignmentColumn);
	}
    }
    
    public Rms rms() {
	if (hasGaps()) throw new RuntimeException("Cannot calculate RMS for AtomAlignment with gaps");
	return new Rms(this);
    }

    public Atom atomAt(int coulumnIndex, int rowIndex) {
	return ((AtomAlignmentColumn) columnAt(coulumnIndex)).atomAt(rowIndex);
    }
}
