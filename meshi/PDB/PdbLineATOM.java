package meshi.PDB;

public class PdbLineATOM extends PdbLineFilter {
    public PdbLineATOM () {}
    public boolean  acceptPdbLine(PdbLine line) {
	return line.isAnAtom();
    }
}
