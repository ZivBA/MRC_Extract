package meshi.util.file;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import meshi.util.CommandList;
import meshi.util.Key;

public class MeshiWriter extends PrintWriter {
    private String path = "unknown";
    public MeshiWriter(String fileName) throws IOException{
	    super(new BufferedWriter(new FileWriter(fileName)));
    }
    
    public MeshiWriter(CommandList commands, Key patheKey, Key nameKey) throws Exception {
	super(getFile(commands,patheKey,nameKey));
	path = (getFile(commands,patheKey,nameKey)).getAbsolutePath() ;
    }

    public MeshiWriter(CommandList commands, Key patheKey, Key nameKey, String postFix) throws Exception {
	super(getFile(commands,patheKey,nameKey, postFix));
	path = (getFile(commands,patheKey,nameKey, postFix)).getAbsolutePath() ;
    }


    private static File getFile(CommandList commands, Key patheKey, Key nameKey) {
	return getFile(commands,patheKey,nameKey, "");
    }
    private static File getFile(CommandList commands, Key patheKey, Key nameKey, String postFix) {
	return new File(commands.firstWord(patheKey).secondWord(),
			commands.firstWord(nameKey).secondWord()+postFix);
    }

    public String toString() {
	return path;
    }
}

