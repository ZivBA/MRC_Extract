package meshi.applications.prediction.analysis;
import java.io.File;
import java.util.Iterator;

public class ExpAlignment extends ModelDataList {
    public ExpAlignment(File directory) {
	super("results from "+directory);
	if (! directory.isDirectory()) throw new RuntimeException(directory+" not a directory");
	String[] fileNames = directory.list();
	for (int i = 0; i < fileNames.length; i++) {
	    String name = fileNames[i];
	    if (name.startsWith("result") & name.endsWith(".pdb")) {
		Model temp = new Model(new File(directory,name));
		if (temp.valid())
		    add(new ModelData(temp));
	    }
	}
    }
    
    public void tabulate(String[] keys) {
	for (Iterator modelsData = iterator(); modelsData.hasNext();) {
	    ModelData modelData = (ModelData) modelsData.next();
	    modelData.tabulate(keys);
	}
    }
	
    public String getXY(String keyX, String keyY) {
	return getXY(keyX, keyY, 0, false);
    }

    public String getXY(String keyX, String keyY, double yShift, boolean perResidue) {
	String out = "";
	for (Iterator modelsData = iterator(); modelsData.hasNext();) {
	    ModelData modelData = (ModelData) modelsData.next();
	    out += modelData.getXY(keyX, keyY,yShift,perResidue);
	}
	return out;
    }
 
    public String GdtVsSolvate() {
	return GdtVsSolvate(0);
    }

    public String GdtVsSolvate(double yShift) {
	String out = "";
	for (Iterator modelsData = iterator(); modelsData.hasNext();) {
	    ModelData modelData = (ModelData) modelsData.next();
	    out += modelData.GdtVsSolvate(yShift);
	}
	return out;
    }

    public double getMinData(String key, boolean perResidue) {
	double min = 10000000;
	for (Iterator models = iterator();models.hasNext();) {
	    ModelData model = (ModelData) models.next();
	    double temp = model.getDoubleValue(key, perResidue);
	    if (temp < min) min = temp;
	}
	return min;
    }

    public String toString() {
	String out =  "ExpAlignment "+comment()+":\n";
	for (Iterator modelsData = iterator();modelsData.hasNext();)
	    out += "\t"+modelsData.next().toString()+"\n";
	return out;
    }
}
    

