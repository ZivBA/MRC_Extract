package meshi.molecularElements.hydrogenBonds;

import meshi.energy.solvate.SolvateParametersList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.AtomList;

/**
 * A subtype of DahiyatHydrogenBondList, that is using a more relaxed angle parameters. 
 */

public class HydrogenBondDahiyatForScmodList extends HydrogenBondDahiyatList {
	
    public HydrogenBondDahiyatForScmodList(DistanceMatrix dm, AtomList atomList,SolvateParametersList parameters) {
    	super(dm, atomList, parameters);
    	angleParameters = new DahiyatLowAccuracyAngleParameters();  
    }

}
