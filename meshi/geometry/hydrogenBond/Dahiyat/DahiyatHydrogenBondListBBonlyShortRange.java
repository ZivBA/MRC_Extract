package meshi.geometry.hydrogenBond.Dahiyat;

import meshi.geometry.DistanceMatrix;
import meshi.geometry.hydrogenBond.AbstractHydrogenBond;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;

public class DahiyatHydrogenBondListBBonlyShortRange extends DahiyatHydrogenBondList {
	
	public final static int SHORT_RANGE_RES_DIS = 6;
	
	public DahiyatHydrogenBondListBBonlyShortRange() {}

	public DahiyatHydrogenBondListBBonlyShortRange(DistanceMatrix dm, AtomList atomList,
			DahiyatParametersInterface parameters) {
		super(dm, atomList, parameters);
	}

    protected AbstractHydrogenBond createHBfromPolars(Atom atom1,Atom atom2) {
    	if ((atom1.isBackbone) && (atom2.isBackbone) && 
    			((atom1.residueNumber()-atom2.residueNumber())<SHORT_RANGE_RES_DIS) &&
    			((atom1.residueNumber()-atom2.residueNumber())>-SHORT_RANGE_RES_DIS))
    		return super.createHBfromPolars(atom1, atom2);    		
    	else
    		return null;
    }
	
}