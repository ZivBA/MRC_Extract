package meshi.applications.prediction.analysis;
import java.io.File;
import java.util.Iterator;
import java.util.StringTokenizer;

import meshi.util.MeshiList;
import meshi.util.filters.Filter;

public class ExpTarget extends MeshiList {
    public ExpTarget(File directory) {
	super(new isExpAlignment());
	setComment("results from "+directory);
	if (! directory.isDirectory()) throw new RuntimeException(directory+" not a directory");
	String[] fileNames = directory.list();
	for (int i = 0; i < fileNames.length; i++) {
	    String name = fileNames[i];
	    if (name.startsWith("0")) {
		ExpAlignment temp = new ExpAlignment(new File(directory,name));
		if (temp.size() > 0)
		    add((Object) temp);
	    }
	}
    }
    
    public void tabulate(String[] keys) {
	System.out.println("TARGET "+comment());
	ExpAlignment alignment = (ExpAlignment) elementAt(0);
	ModelData modelData =  (ModelData) alignment.elementAt(0);
	System.out.printf("%-10s ","length");
	for (int i = 0; i < keys.length; i++) { 
	    String key = keys[i];
	    System.out.printf("%-9s ",key.substring(3));
	}
	System.out.println("path");

	for (Iterator alignments = iterator();alignments.hasNext();) {
	    alignment = (ExpAlignment) alignments.next();
	    alignment.tabulate(keys);
	}
    }
	
    public String getXY(String keyX, String keyY, double yShift,boolean perResidue) {
	String out = "";
	if (yShift == 0) {
	    for (Iterator alignments = iterator();alignments.hasNext();) {
		ExpAlignment alignment = (ExpAlignment) alignments.next();
		out += alignment.getXY(keyX,keyY,yShift,perResidue);
	    }
	    return out;
	}
	else {
	    double minY = 10000000;
	    for (Iterator alignments = iterator();alignments.hasNext();) {
		ExpAlignment alignment = (ExpAlignment) alignments.next();
		double temp = alignment.getMinData(keyY,perResidue);
		if (temp < minY) minY = temp;
	    }
	    for (Iterator alignments = iterator();alignments.hasNext();) {
		ExpAlignment alignment = (ExpAlignment) alignments.next();
		out += alignment.getXY(keyX,keyY,yShift- minY,perResidue);
	    }
	    return out;
	}
    }
	
    public String toString() {
	String out =  "ExpTarget "+comment()+"\n";
	for (Iterator alignments = iterator();alignments.hasNext();) {
	    StringTokenizer lines = new StringTokenizer(alignments.next().toString(),"\n");
	    while (lines.hasMoreTokens())
		out += "\t"+lines.nextToken()+"\n";
	}
	return out;
    }

    public static class isExpAlignment implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof ExpAlignment);
	}
    }
 }
