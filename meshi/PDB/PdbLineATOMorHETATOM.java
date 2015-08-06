package meshi.PDB;

public class PdbLineATOMorHETATOM extends PdbLineFilter {
    public PdbLineATOMorHETATOM () {}
    public boolean  acceptPdbLine(PdbLine line) {
	return line.isAnAtomOrHeteroAtom();
    }
}
