package meshi.sequences;
import java.util.Iterator;

import meshi.molecularElements.DummyResidue;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.util.filters.Filter;

public class ResidueAlignment extends Alignment{
    /**
     * Empty alignment.
     **/
    public ResidueAlignment() {
	super(new IsAlignmentColumn());
    }

    /**
     * A trevial alignment of two protein objects which are assumed to be two models of the same protein.
     **/
    public ResidueAlignment(Protein protein1, Protein protein2) {
	this();
	comments.add(protein1.name());
	comments.add(protein2.name());
	Iterator residues1 = protein1.residues().iterator();
	Iterator residues2 = protein2.residues().iterator();
	while (residues1.hasNext() & residues2.hasNext()) {
	    Residue residue1 = (Residue) residues1.next();
	    Residue residue2 = (Residue) residues2.next();
	    ResidueAlignmentColumn column = new ResidueAlignmentColumn(residue1,residue2);
	    add(column);
	}
	if (residues1.hasNext())
	    while (residues1.hasNext()) {
		Residue residue1 = (Residue) residues1.next();
		add(new ResidueAlignmentColumn(new ResidueAlignmentCell(residue1),
					       new ResidueAlignmentCell(new DummyResidue(residue1.number))));
	    }
	if (residues2.hasNext())
	    while (residues2.hasNext()) {
		Residue residue2 = (Residue) residues2.next();
		add(new ResidueAlignmentColumn(new ResidueAlignmentCell(new DummyResidue(residue2.number)),
					       new ResidueAlignmentCell(residue2)));
	    }
    }


    public String toString() {
	return (new SequenceAlignment(this)).toString();
    }
       
    private static class IsAlignmentColumn implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof ResidueAlignmentColumn);
	}
    }
}
