package meshi.applications.minimizeProtein;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
/**
 * An extended atoms model of a protein.
 * This model includes all heavy atoms and most hydrogen atoms that take part in hydrogen bonds.
 **/ 
public class ExtendedAtomsProtein extends Protein {
    /**
     * The parameter is expected to be a PDB formatted file.
     **/ 
    public ExtendedAtomsProtein(String fileName) {
	super((new AtomList(fileName)), new ResidueExtendedAtoms());
    }

    public ExtendedAtomsProtein(String fileName,int addAtomsFlag) {
	super((new AtomList(fileName)), new ResidueExtendedAtoms(addAtomsFlag));
    }

    public ExtendedAtomsProtein(AtomList atoms,int addAtomsFlag) {
	super(atoms, new ResidueExtendedAtoms(addAtomsFlag));
    }

}
