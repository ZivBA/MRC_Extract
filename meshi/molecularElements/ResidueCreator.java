package meshi.molecularElements;
public interface ResidueCreator {
    public Residue create(AtomList atoms, int residueNumber, int mode);
    public Residue create(String name, int residueNumber, int mode, double x, double y, double z);
}
