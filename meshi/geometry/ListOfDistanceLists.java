package meshi.geometry;

import meshi.util.MeshiList;
import meshi.util.filters.Filter;

public class ListOfDistanceLists extends MeshiList {

    protected ListOfDistanceLists() {
           super(new IsDistanceList(), 2);
       }

 public static class IsDistanceList implements Filter {
    public boolean accept(Object obj) {
        return (obj instanceof DistanceList);
        }
    }

    public final DistanceList distanceListAt(int i) { return (DistanceList) elementAt(i);}
}

