package meshi.PDB;

public class PdbLineSEQRES extends PdbLineFilter {
    public boolean  acceptPdbLine(PdbLine line) {
	return line.isSEQRES();
    }
}
