package meshi.energy.hydrogenBondsPlane;

import meshi.geometry.DistanceList;

public class CNdistanceList extends DistanceList  {
    protected CNdistanceList(int capacity) {
           super(capacity);
       }
    protected CNdistanceList() {
           super();
       }

    private IsGoodCNPair isGoodCNPair = new IsGoodCNPair();

    public boolean selectionAdd(Object object){
     return super.selectionAdd(object,isGoodCNPair);
    }

}

