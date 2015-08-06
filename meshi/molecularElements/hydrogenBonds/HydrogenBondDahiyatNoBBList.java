package meshi.molecularElements.hydrogenBonds;

import meshi.energy.solvate.SolvateParametersList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;

/**
 * A subtype of DahiyatHydrogenBondList, that is keeping any H-bonds that are NOT between two backbone atoms, i.e. backbone-sidechain and sidechain-sidechain 
 **/

public class HydrogenBondDahiyatNoBBList extends HydrogenBondDahiyatList {
    
    public HydrogenBondDahiyatNoBBList(DistanceMatrix dm, AtomList atomList,SolvateParametersList parameters) {
	super(dm, atomList, parameters);
    }
    
    protected AbstractHydrogenBond createHBfromPolars(Atom atom1,Atom atom2) {
	if (!atom1.isBackbone || !atom2.isBackbone)
	    return super.createHBfromPolars(atom1, atom2);    
	else
	    return null;
    }

}
