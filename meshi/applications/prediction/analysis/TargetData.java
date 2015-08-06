package meshi.applications.prediction.analysis;
import java.io.File;
import java.util.StringTokenizer;

import meshi.util.filters.Filter;
import meshi.util.string.StringList;

public class TargetData extends ModelDataList {
    private boolean verbose = false;
    public final Filter fileFilter;

    public TargetData(File directory, int startSimulation, int endSimulation ) {
	super("Analysis of "+directory.getPath());
	System.out.print("analysing "+directory+"\t");
	fileFilter = new FileFilter(startSimulation, endSimulation);
	File[] alignments = directory.listFiles();
	for (int i = 0; i < alignments.length; i++) {
	    File alignment = alignments[i];
	    if (alignment.isDirectory()) {
               if (verbose) {
		       System.out.print("analysing "+alignment+"  ");
		       if (i%5 == 0) System.out.println();
	       }
	       else { 
			System.out.print(".");
			if (i%70 == 0) System.out.println();
	       }
 	       File[] files = alignment.listFiles();
	       for (int j = 0; j < files.length; j++) {	
		       File file = files[j];
	    	       if(fileFilter.accept(file)) {
		       		if (verbose) System.out.println("adding "+file);
		       		add((new Model(file)).modelData);
		       }
	       }
	    }
	}
	System.out.println();
    }

    private static class FileFilter implements Filter {
	public final int high, low;

	public FileFilter(int low, int high) {
	    this.high = high;
	    this.low = low;
	}

	public boolean accept(Object obj) {
	    File file = (File) obj;
	    StringList fileName = new StringList(new StringTokenizer(file.getName(),"."));
	    if (fileName.size() != 3) return false; 
	    if (! fileName.stringAt(0).equals("result")) return false;
	    int number = Integer.valueOf(fileName.stringAt(1));
	    if ((number <low) | (number > high)) return false;
	    return true;
	}
    }

    private static class IsModel implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof Model);
	}
    }

}	    
