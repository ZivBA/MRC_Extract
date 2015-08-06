package meshi.PDB;
import java.io.File;

import meshi.util.filters.Filter;
public class isAPdbFile implements Filter {
    public boolean accept(Object obj) {
	if (!(obj instanceof File)) 
	    throw new RuntimeException("weird input to isAPdbFile.accept\n"+
				       obj);
	File file = (File) obj;
	String name = file.getName();
	if (name.endsWith(".pdb")) return true;
	if (name.endsWith(".pdb.gz")) return true;
 	if (name.endsWith(".ent.gz")) return true;
 	if (name.endsWith(".ent")) return true;
	return false;
    }
}
