package meshi.energy.hydrogenBondsPlane;

import meshi.util.MeshiList;

/**
 **/
public class CNtwoDistancesList extends MeshiList  {
    protected CNtwoDistancesList(int capacity) {
           super(new CNtwoDistances.IsCNtwoDistances(),capacity);
       }

    protected CNtwoDistancesList() {
           super(new CNtwoDistances.IsCNtwoDistances());
       }
}

