package meshi.energy.torsionSpaceMinimization;

import java.util.Vector;

import meshi.geometry.Angle;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionList;
import meshi.molecularElements.Atom;

/**
 * This class encapsulates the No. 3 atom in each torsion, with all the data that is relevant to it.
 * It is supposed to be created and run from "TorsionMinimizationTree" 
 **/
public class TorsionMinimizationElement {
	private DistanceMatrix dm;
	private final Atom a1,a2,a3;
	private final Atom mainTorsionAtom;
	private final double angle;
	private final double bond;
	private final boolean isElementAroundSignificantTorsion; // i.e. around any IUPAC names torsion that is not OMEGA
	private double[] torsionCoors;      // 0 - torsion value    1 - torsion gradient
	private double[] grad = new double[6];
	private TorsionMinimizationElement sonMain = null;
	private Distance disMain;
	private Distance axis;
	
	// All the data that is needed to know about the other branches
	private class DataOtherBranch {
		Atom otherAtom;
		double angle;
		double bond;
		double torsionDiffFromMain;
		Distance dis = null;
		TorsionMinimizationElement son = null;
	}
	private Vector<DataOtherBranch> otherTorsionAtoms;
	
	protected TorsionMinimizationElement(Atom a1, Atom a2, Atom a3, Atom mainTorsionAtom, DistanceMatrix dm) {
		this.dm = dm;
		this.a1 = a1;
		this.a2 = a2;
		this.a3 = a3;
		this.mainTorsionAtom = mainTorsionAtom;
		otherTorsionAtoms = new Vector<DataOtherBranch>();
		torsionCoors = new double[2];
		axis = dm.distance(a2,a3);
		Angle angleInstance = new Angle(a2,a3,mainTorsionAtom,dm,false);
		angle = angleInstance.angle();
		disMain = dm.distance(a3, mainTorsionAtom);
		bond = disMain.distance();
		Torsion torsion = new Torsion(new Angle(a1,a2,a3,dm,false), angleInstance, dm);
		torsionCoors[0] = torsion.torsion();
		if (TorsionList.isNameInIUPAC(torsion) && !TorsionList.isOmega(torsion))
			isElementAroundSignificantTorsion = true;
		else
			isElementAroundSignificantTorsion = false;
	}

	protected void addBranch(Atom otherAtom, TorsionMinimizationElement tme) {
		DataOtherBranch dataOtherBranch = new DataOtherBranch();
		dataOtherBranch.otherAtom = otherAtom;
		Angle angleInstance = new Angle(a2,a3,otherAtom,dm,false);
		dataOtherBranch.angle = angleInstance.angle();
		dataOtherBranch.dis = dm.distance(a3, otherAtom);
		dataOtherBranch.bond = dataOtherBranch.dis.distance();
		Torsion torsion = new Torsion(new Angle(a1,a2,a3,dm,false), angleInstance, dm);
		dataOtherBranch.torsionDiffFromMain = torsion.torsion() - torsionCoors[0];
		dataOtherBranch.son = tme;
		otherTorsionAtoms.add(dataOtherBranch);
	}
	
	protected void setSonMain(TorsionMinimizationElement tme) {
		sonMain = tme;
	}
	
	protected void buildTorsionInProt(TempCalcForthAtom tempCalcForthAtom) {
		tempCalcForthAtom.setForthAtom(bond, angle, torsionCoors[0], a1, a2, a3, mainTorsionAtom);
		for (DataOtherBranch dob : otherTorsionAtoms)
			tempCalcForthAtom.setForthAtom(dob.bond, dob.angle, torsionCoors[0] + dob.torsionDiffFromMain, a1, a2, a3, dob.otherAtom);
	}
	
	protected void calcDericatives() {
		//mainTorsion.update();
		//for (DataOtherBranch dob : otherTorsionAtoms)
		//	dob.torsion.update();
		grad[0] = grad[1] = grad[2] = grad[3] = grad[4] = grad[5] = 0.0;
		// The gradient from the main branch.
		if (sonMain==null) {  // The main branch is a terminal
			grad[0]+=mainTorsionAtom.fx();
			grad[1]+=mainTorsionAtom.fy();
			grad[2]+=mainTorsionAtom.fz();
			grad[3]+=  /* rotation around X */
				mainTorsionAtom.fz()* disMain.dy() - mainTorsionAtom.fy()* disMain.dz();
			grad[4]+=  /* rotation around Y */
				mainTorsionAtom.fx()*disMain.dz() - mainTorsionAtom.fz()* disMain.dx();
			grad[5]+=  /* rotation around Z */
				mainTorsionAtom.fy()* disMain.dx() - mainTorsionAtom.fx()*disMain.dy();
		} 
		else {  // There is a continuation from the main branch
			grad[0] += mainTorsionAtom.fx() + sonMain.grad[0];
			grad[1] += mainTorsionAtom.fy() + sonMain.grad[1];
			grad[2] += mainTorsionAtom.fz() + sonMain.grad[2];
			grad[3] += disMain.dy()*(mainTorsionAtom.fz() + sonMain.grad[2]) - disMain.dz()*(mainTorsionAtom.fy() + sonMain.grad[1]) +
				sonMain.grad[3];
			grad[4] += disMain.dz()*(mainTorsionAtom.fx() + sonMain.grad[0]) - disMain.dx()*(mainTorsionAtom.fz() + sonMain.grad[2]) +
			    sonMain.grad[4];
			grad[5] += disMain.dx()*(mainTorsionAtom.fy() + sonMain.grad[1]) - disMain.dy()*(mainTorsionAtom.fx() + sonMain.grad[0]) +
			    sonMain.grad[5];
		}
		
		// The gradients from the other branches
		for (DataOtherBranch dob : otherTorsionAtoms) {
			if (dob.son==null) {  // The main branch is a terminal
				grad[0]+=dob.otherAtom.fx();
				grad[1]+=dob.otherAtom.fy();
				grad[2]+=dob.otherAtom.fz();
				grad[3]+=  /* rotation around X */
					dob.otherAtom.fz()* dob.dis.dy() - dob.otherAtom.fy()* dob.dis.dz();
				grad[4]+=  /* rotation around Y */
					dob.otherAtom.fx()*dob.dis.dz() - dob.otherAtom.fz()* dob.dis.dx();
				grad[5]+=  /* rotation around Z */
					dob.otherAtom.fy()* dob.dis.dx() - dob.otherAtom.fx()*dob.dis.dy();
			} 
			else {  // There is a continuation from the main branch
				grad[0] += dob.otherAtom.fx() + dob.son.grad[0];
				grad[1] += dob.otherAtom.fy() + dob.son.grad[1];
				grad[2] += dob.otherAtom.fz() + dob.son.grad[2];
				grad[3] += dob.dis.dy()*(dob.otherAtom.fz() + dob.son.grad[2]) - dob.dis.dz()*(dob.otherAtom.fy() + dob.son.grad[1]) +
					dob.son.grad[3];
				grad[4] += dob.dis.dz()*(dob.otherAtom.fx() + dob.son.grad[0]) - dob.dis.dx()*(dob.otherAtom.fz() + dob.son.grad[2]) +
				    dob.son.grad[4];
				grad[5] += dob.dis.dx()*(dob.otherAtom.fy() + dob.son.grad[1]) - dob.dis.dy()*(dob.otherAtom.fx() + dob.son.grad[0]) +
				    dob.son.grad[5];
			}			
		}
		
		// The gradient according to the main torsion
		torsionCoors[1] = grad[3]*axis.dDistanceDx() + grad[4]*axis.dDistanceDy() + grad[5]*axis.dDistanceDz();
	}
	
	public boolean isElementAroundSignificantTorsion() {
		return isElementAroundSignificantTorsion;
	}
	
	public double[] getCoors() {
		return torsionCoors;
	}
	
}
