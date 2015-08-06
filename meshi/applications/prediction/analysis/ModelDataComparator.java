package meshi.applications.prediction.analysis;
import java.util.Comparator;

public class ModelDataComparator implements Comparator {
    private String key;
    public ModelDataComparator(String key) {
	this.key = key;
    }
    public int compare(Object obj1, Object obj2) {
	ModelData model1 = (ModelData) obj1;
	ModelData model2 = (ModelData) obj2;
	Double value1 = new Double((String) model1.data.get(key));
	Double value2 = new Double((String) model2.data.get(key));
	return value1.compareTo(value2);
    }
}
	
