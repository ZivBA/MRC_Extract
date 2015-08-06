package meshi.sequences;
import meshi.molecularElements.Atom;

public class AtomAlignmentCell extends AlignmentCell {
    public AtomAlignmentCell(Atom atom, String comment) {
	super(atom, atom.number(), comment);
    }
    
    public Atom atom() {return (Atom) obj;}
    public boolean gap() {
	return (obj == null);
    }
}
