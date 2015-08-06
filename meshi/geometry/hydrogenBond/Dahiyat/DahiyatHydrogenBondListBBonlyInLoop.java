package meshi.geometry.hydrogenBond.Dahiyat;

import meshi.geometry.DistanceMatrix;
import meshi.geometry.hydrogenBond.AbstractHydrogenBond;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;

public class DahiyatHydrogenBondListBBonlyInLoop extends DahiyatHydrogenBondList {

	private int loopStart=-1;
	private int loopEnd=-1;
	
	
	public DahiyatHydrogenBondListBBonlyInLoop() {}

	public DahiyatHydrogenBondListBBonlyInLoop(DistanceMatrix dm, AtomList atomList,
			DahiyatParametersInterface parameters, int loopStart, int loopEnd) {
		super(dm, atomList, parameters);
		this.loopEnd = loopEnd;
		this.loopStart = loopStart;
	}

    protected AbstractHydrogenBond createHBfromPolars(Atom atom1,Atom atom2) {
    	if ((atom1.isBackbone) && (atom2.isBackbone) && 
    			(atom1.residueNumber()>=loopStart) && (atom1.residueNumber()<=loopEnd) &&
    			(atom2.residueNumber()>=loopStart) && (atom2.residueNumber()<=loopEnd))
    		return super.createHBfromPolars(atom1, atom2);    		
    	else
    		return null;
    }
	
	
}
