package meshi.geometry.hydrogenBond.Dahiyat;

import meshi.geometry.DistanceMatrix;
import meshi.geometry.hydrogenBond.AbstractHydrogenBond;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;

public class DahiyatHydrogenBondListBBonlyLongRange extends DahiyatHydrogenBondList {
		
	public DahiyatHydrogenBondListBBonlyLongRange() {}

	public DahiyatHydrogenBondListBBonlyLongRange(DistanceMatrix dm, AtomList atomList,
			DahiyatParametersInterface parameters) {
		super(dm, atomList, parameters);
	}

    protected AbstractHydrogenBond createHBfromPolars(Atom atom1,Atom atom2) {
    	if ((atom1.isBackbone) && (atom2.isBackbone) && 
    			(((atom1.residueNumber()-atom2.residueNumber())>=DahiyatHydrogenBondListBBonlyShortRange.SHORT_RANGE_RES_DIS) ||
    			((atom1.residueNumber()-atom2.residueNumber())<=-DahiyatHydrogenBondListBBonlyShortRange.SHORT_RANGE_RES_DIS)))
    		return super.createHBfromPolars(atom1, atom2);    		
    	else
    		return null;
    }
	
}
