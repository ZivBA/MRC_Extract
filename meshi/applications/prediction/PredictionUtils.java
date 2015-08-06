package meshi.applications.prediction;
import meshi.util.CommandList;
import meshi.util.Key;
import meshi.util.MeshiProgram;
import meshi.util.file.MeshiWriter;

/**
 * Where useful static methods are found. 
 * This is NOT the right place for highly sofisticate algorithms. 
 **/
public class PredictionUtils {
    /**
     * Opens a file for writing (warped by MeshiWriter object) with a name and location specified
     * by the commands and a .<seed>.pdb suffix.
     **/
    public static MeshiWriter newPdbWriter(CommandList commands, Key pathKey, Key nameKey) {
	try {
	    String suffix = "."+MeshiProgram.seed()+".pdb";
	    return new MeshiWriter(commands, pathKey, nameKey, suffix);
	}
	catch(Exception ex) {
	    throw new RuntimeException("A problem in opening output file with key words "+
				       pathKey+" and "+nameKey+"\n"+ex);
	}
    }
}
