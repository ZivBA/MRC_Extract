package meshi.energy.hydrogenBond;

import meshi.geometry.DistanceList;
import meshi.util.filters.Filter;
/**
 **/
public class HBdistanceList extends DistanceList {
    protected HBdistanceList(int capacity) {
           super(capacity);
       }
    protected HBdistanceList() {
           super();
       }

    public HBdistanceList(Filter fi) {
        super(fi);
    }

    private GoodResiduesForHB goodResiduesForHB = new GoodResiduesForHB();

    public boolean selectionAdd(Object object){
     return this.selectionAdd(object,goodResiduesForHB);
    }
    
}

