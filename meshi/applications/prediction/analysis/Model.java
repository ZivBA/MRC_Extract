package meshi.applications.prediction.analysis;
import java.io.File;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.StringTokenizer;

import meshi.PDB.PdbReader;
import meshi.molecularElements.Protein;
import meshi.util.string.StringList;

public class Model {
    public final Hashtable dictionary, data;
    public final ModelData modelData;
    private static DictionaryLineFilter dictionaryLineFilter = new DictionaryLineFilter();
    public final File file;
    public final Protein protein;
    public final Exception exception;
    public final int numberOfResidues;
    public Model(String fileName) {
	    this(new File(fileName));
    }

    public Model(File file) {
	Exception ex1 = null;
	Protein protein1 = null;
	try {
	    protein1 = new Protein(new PdbReader(file));
	}
	catch (Exception ex) {ex1 = ex;}
	if (ex1 == null) {
	    dictionary = buildDictionary(file);
	    data = buildData(file, dictionary);
	    exception = null;
	    numberOfResidues = protein1.residues().numberOfNonDummyResidues();
	}
	else {
	    exception = ex1;
	    dictionary = null;
	    data = null;
	    numberOfResidues = -1;
	}
	this.file = file;
	modelData = new ModelData(this);
	protein = protein1;
    }

    public String toString() {
	if (!valid()) return "Invalid Model";
	return "model of "+protein;
    }

    public boolean valid() {
	return (exception == null);
    }

    private static Hashtable buildDictionary(File file) {
	Hashtable out = new Hashtable();
	StringList tempList = new StringList(file);
	StringList dictionaryLines = tempList.filter(dictionaryLineFilter);
	for (Iterator lines = dictionaryLines.iterator(); lines.hasNext();){
	    String line = (String) lines.next();
	    StringList words = new StringList(new StringTokenizer(line));
	    String key = words.stringAt(2);
	    String value = line.substring(line.indexOf(key));
	    out.put(key,value);
	}
	return out;
    }

    private static Hashtable buildData(File file,Hashtable dictionary) {
	KeyLineFilter[] keyFilters = { new KeyLineFilter("T1"),
				       new KeyLineFilter("T2"),
				       new KeyLineFilter("T3"),
				       new KeyLineFilter("T4") };
	KeyLineFilter[] valueFilters = { new KeyLineFilter("V1"),
					 new KeyLineFilter("V2"),
					 new KeyLineFilter("V3"),
					 new KeyLineFilter("V4") };
	Hashtable out = new Hashtable();
	StringList tempList = new StringList(file);

	for (int i = 0; i < keyFilters.length;i++) {
	    KeyLineFilter kf = keyFilters[i];
	    KeyLineFilter vf = valueFilters[i];
	    String keysLine = null;
	    String valuesLine = null;
	    for (Iterator lines = tempList.iterator(); lines.hasNext();) {
		String line = (String) lines.next();
		if (kf.accept(line)) keysLine = line;
		if (vf.accept(line)) valuesLine = line;
		if ((keysLine != null) & (valuesLine != null)) break;
	    }
	    if ((keysLine == null)) throw new RuntimeException("Did not find key "+kf.key+" in "+file); 
	    if ((valuesLine == null)) throw new RuntimeException("Did not find key "+vf.key+" in "+file);
	    
	    StringTokenizer keys = new StringTokenizer(keysLine);
	    StringTokenizer values = new StringTokenizer(valuesLine);
	    while (keys.hasMoreTokens() & values.hasMoreTokens()) {
		String key = kf.key+"."+keys.nextToken();
		String value = values.nextToken();
		    out.put(key, value);
	    }
	}
	return out;
    }    

}	 

	
	
	
