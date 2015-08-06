package meshi.applications.prediction.analysis;
import java.util.Hashtable;
import java.util.Iterator;

import meshi.util.SortableMeshiList;
import meshi.util.filters.Filter;

public class ModelDataList extends SortableMeshiList {
    private boolean verbose = true;
    public double enrichmentRms = -1;

    public ModelDataList(String comment){
	super(new IsModelData());
	setComment(comment);
    }
    
    public Hashtable calculateAverages(String[] keys) {
	if (verbose) System.out.println("calculating averages ");
	Hashtable out = new Hashtable();
	double[] avg = new double[keys.length];
	double[] avg2 = new double[keys.length];
	for (int i = 0; i < keys.length; i++) 
	    avg[i] = avg2[i] = 0;

	int numberOfValidModels = 0;
	for (Iterator models = iterator(); models.hasNext();) {
	    ModelData modelData = (ModelData) models.next();
	    if (modelData.valid) {
		numberOfValidModels++;
		if (verbose) System.out.print("adding "+modelData.file+"  ");
		for (int iKey = 0; iKey < keys.length; iKey++) {
		    String key = keys[iKey];
		    double doubleValue;
		    try {
			doubleValue = Double.parseDouble((String) modelData.data.get(key));
		    }
		    catch (Exception ex) { 
			throw new RuntimeException("Failed to parse "+modelData.data.get(key)+"\n"+
						   "key = "+key+"\n"+ex+"\n"+
						   modelData.file.getPath());
		    }
		    if (verbose) System.out.print(key+"="+doubleValue+"; ");
			
		    avg[iKey] += doubleValue;
		    avg2[iKey] += doubleValue*doubleValue;

		}
	    }
	    if (verbose) System.out.println();
	}
	for (int iKey = 0; iKey < keys.length; iKey++) {
	    avg[iKey] /= numberOfValidModels;
	    avg2[iKey] /= numberOfValidModels;
	}
	    
	for (int iKey = 0; iKey < keys.length; iKey++) {
	    ListData listData = new ListData(avg[iKey], 
			    		     Math.sqrt(avg2[iKey]-avg[iKey]*avg[iKey]),
					     size(),
					     keys[iKey]+" - "+comment());
	    out.put(keys[iKey],listData);
	}
	return out;
    }

    public ModelDataList extractLowest(String key, int n) {
	    ModelDataList tempList = (ModelDataList) filter(new ValidFilter(), new ModelDataList(comment()));
	    ModelDataList lowestList = null;
	    try {
		    lowestList = (ModelDataList) tempList.extractLowest(new ModelDataComparator(key), n, 
			                     new ModelDataList("Lowest "+key+" of "+comment()));
	    }
	    catch (Exception ex) { System.out.println(ex);}
	    if (lowestList == null) return new ModelDataList(comment());
	    return lowestList;

    }

    public double lowestRms() {
	    double out = 10000;
	    for (Iterator models = iterator(); models.hasNext();) {
		    ModelData modelData = (ModelData) models.next();
		    if (modelData.valid) {
		   	double rms = Double.valueOf((String)modelData.data.get("rms1")); 
		        if (out > Double.valueOf(rms)) out = rms;
		    }
	    }
	    return out;
    }

    private static class ValidFilter implements Filter {
	      public boolean accept(Object obj) {
		      ModelData md = (ModelData) obj;
		      return md.valid;
	      }
    }
    private static class IsModelData implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof ModelData);
	}
    }

} 
