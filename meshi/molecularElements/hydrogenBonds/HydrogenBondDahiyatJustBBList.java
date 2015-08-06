package meshi.molecularElements.hydrogenBonds;

import meshi.energy.solvate.SolvateParametersList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;

/**
 * A subtype of DahiyatHydrogenBondList, that is keeping only backbone hydrogen bonds (i.e. backbone to backbone) 
 */

public class HydrogenBondDahiyatJustBBList extends HydrogenBondDahiyatList {
	
    public HydrogenBondDahiyatJustBBList(DistanceMatrix dm, AtomList atomList,SolvateParametersList parameters) {
    	super(dm, atomList, parameters);
    }
    
    protected AbstractHydrogenBond createHBfromPolars(Atom atom1,Atom atom2) {
    	if ((atom1.isBackbone) && (atom2.isBackbone))
    		return super.createHBfromPolars(atom1, atom2);    		
    	else
    		return null;
    }

}
