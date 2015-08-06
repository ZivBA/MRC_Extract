package meshi.geometry;

import meshi.util.SortableMeshiList;
import meshi.util.filters.Filter;

public class DistanceList extends SortableMeshiList {

    protected DistanceList(int capacity) {
           super(new IsDistance(), capacity);
       }

    public DistanceList() {
        super(new IsDistance());
    }

    public DistanceList(Filter fi) {
        super(fi);
    }
    public final Distance distanceAt(int i) { return (Distance) elementAt(i);}

    private static class IsDistance implements Filter {
    public boolean accept(Object obj) {
        return (obj instanceof Distance);
        }
    }

    private IsDistance filter = new IsDistance();

    public boolean selectionAdd(Object object){
      return this.selectionAdd(object,filter);
     }

     public boolean selectionAdd(Object object, Filter filter){
        if (filter.accept(object)) {
            fastAdd(object);
            return true;
        }
         else return false;
     }

    public Object[] internalArray() {return internalArray;}
}

