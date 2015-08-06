 package meshi.PDB;
import meshi.util.filters.Filter;
public abstract class PdbLineFilter implements Filter {
    public boolean accept(Object obj) {
	return acceptPdbLine((PdbLine) obj);
    }
    public abstract boolean acceptPdbLine(PdbLine line);
}
