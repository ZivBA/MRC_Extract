
package meshi.molecularElements;
public class DummyResidue extends Residue {
    public DummyResidue(int residueNumber) {
	super("UKN",UNK,residueNumber,new AtomList());
	dummy = true;
    }

    public DummyResidue(DummyResidue old) {
	this(old.number);	
    }

    public Atom ca() { return null;}
}
