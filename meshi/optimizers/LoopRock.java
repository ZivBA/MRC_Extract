package meshi.optimizers;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.evRot1BB.EVRot1BBCreator;
import meshi.energy.evRot1BB.EVRot1BBEnergy;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBondsAngle.HBondsPunishOHNAngleCreator;
import meshi.energy.hydrogenBondsAngle.HbondsPunishHOCAngleCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsCreator;
import meshi.energy.hydrogenBondsPlane.HydrogenBondsPlaneCreator;
import meshi.energy.linearRG.LinearRgCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.softExcludedVol.SoftExcludedVol;
import meshi.energy.softExcludedVol.SoftExcludedVolCreator;
import meshi.energy.solvateRot1.SolvateRot1Creator;
import meshi.energy.solvateRot1.SolvateRot1Energy;
import meshi.energy.torsionVal.TorsionValCreator;
import meshi.energy.twoTorsions.TwoTorsionsCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.Helixer;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.rotamericTools.RotamericTools;


public class LoopRock {


private Protein prot;
private int resBegin , resEnd;
private CommandList commands;
private DistanceMatrix dm;
private DunbrackLib lib;
private double[][] pp; 
private double weightEV;


public LoopRock(Protein prot , int resBegin , int resEnd , int nIter , CommandList commands) {
	this(prot , resBegin , resEnd , nIter , 1 , commands);
}


public LoopRock(Protein prot , int resBegin , int resEnd , int nIter , double weightEV, CommandList commands) {
	this.prot = prot;
	this.resBegin = resBegin;
	this.resEnd = resEnd;
	this.weightEV = weightEV;
	this.commands = commands;
	lib = new DunbrackLib(commands, 1.0 , 2);
	rock(nIter);
}

private void rock(int nIter) {
	double angJump = 10.0;  // must divide 360 exactly
	EnergyCreator[] energyCreators = new EnergyCreator[15];
	TotalEnergy energy;
	Minimizer minimizer;
	AbstractEnergy realev;
	AbstractEnergy ev;
	AbstractEnergy sol;
	rot1andFreeze();
	int defrostedCount = 0;
	for (int c=0 ; c<prot.atoms().size() ; c++)
		if (!prot.atoms().atomAt(c).frozen())
			defrostedCount++;
	System.out.println("Number of defrosted atoms: " + defrostedCount);
	double[][] bestPos = new double[defrostedCount][3];
	double bestMinEne = 1e10;

	for ( ; nIter>0 ; nIter--) {
	System.out.println("\nIterations to go: " + nIter + "\n");
	rot1andFreeze();
	energyCreators[0] = new BondCreator();
	energyCreators[1] = new AngleCreator();
	energyCreators[2] = new PlaneCreator();
	energyCreators[3] = new OutOfPlaneCreator();
	energyCreators[4] = new SoftExcludedVolCreator(100.0,2);
	energyCreators[5] = new EVRot1BBCreator(weightEV,pp,1000);
	energyCreators[6] = new SolvateRot1Creator(1,pp,1000);
	energyCreators[7] = new TwoTorsionsCreator(30.0);
	energyCreators[8] = new LinearRgCreator(10.0);
	energyCreators[9] = new TorsionValCreator(1000.0);
	energyCreators[10] = new HydrogenBondsCreator(0.05);
	energyCreators[11] = new HydrogenBondsPairsCreator(0.05,(HydrogenBondsCreator) energyCreators[10]);
	energyCreators[12] = new HBondsPunishOHNAngleCreator(5.0,(HydrogenBondsCreator) energyCreators[10]);
	energyCreators[13] = new HbondsPunishHOCAngleCreator(5.0,(HydrogenBondsCreator) energyCreators[10]);
	energyCreators[14] = new HydrogenBondsPlaneCreator(0.05);
		
	energy = new TotalEnergy(prot, dm, energyCreators, commands);
	realev = energy.getEnergyTerm(new SoftExcludedVol());
	ev = energy.getEnergyTerm(new EVRot1BBEnergy());
	sol = energy.getEnergyTerm(new SolvateRot1Energy());
	
	for (int c=resBegin ; c<= resEnd ; c++) {
		double bestEne = 1e10;
		double bestInd = 36; 
		for (int d=0 ; d< Math.round(360/angJump) ; d++) {
			energy.resetAtomEnergies();
			try {
				energy.update();
			}
			catch (Exception e) {
				throw  new RuntimeException(e);
			}
			ev.evaluateAtoms();
			sol.evaluateAtoms();
			double ene = prot.residue(c).ca().energy();
			if (ene<bestEne) {
				bestEne = ene;
				bestInd = d;
			}
			rotateRes(c,angJump*Math.PI/180.0);			
		}
		System.out.println(bestInd + " " + bestEne);
		rotateRes(c,bestInd*angJump*Math.PI/180.0);		
	}
	// Now Minimizing
	minimizer = new LBFGS(energy, 0.05 , 10000 , 100); 
	try {
		System.out.println(minimizer.minimize());
	} 
	catch (Exception e) {
		throw  new RuntimeException(e);
	}
	if (bestMinEne > ev.evaluate() + sol.evaluate()) {
		bestMinEne = ev.evaluate() + sol.evaluate();
		System.out.println("New best min energy: " + bestMinEne);
		defrostedCount = 0;
		for (int c=0 ; c<prot.atoms().size() ; c++)
			if (!prot.atoms().atomAt(c).frozen()) {
				bestPos[defrostedCount][0] = prot.atoms().atomAt(c).x();
				bestPos[defrostedCount][1] = prot.atoms().atomAt(c).y();
				bestPos[defrostedCount][2] = prot.atoms().atomAt(c).z();
				defrostedCount++;
			}
	}
	}  // Of nIter
	defrostedCount = 0;
	for (int c=0 ; c<prot.atoms().size() ; c++)
		if (!prot.atoms().atomAt(c).frozen()) {
			prot.atoms().atomAt(c).setXYZ(bestPos[defrostedCount][0],bestPos[defrostedCount][1],bestPos[defrostedCount][2]);
			defrostedCount++;
		}
}

private void rot1andFreeze() {
	prot.defrost();
	dm = new DistanceMatrix(prot.atoms(),  5.5, 2.0, 4); 
	pp = RotamericTools.putIntoRot1(prot, dm, lib);
	prot.freeze();
	for (int c=resBegin ; c<= resEnd ; c++) 
		prot.residue(c).atoms().defrost();
	dm = new DistanceMatrix(prot.atoms(),  5.5, 2.0, 4); 	
}

private void rotateRes(int resNum , double ang) {
	AtomList list = prot.residue(resNum).atoms();
	Atom cAtom = find(list,"C");
	Atom nAtom = find(list,"N");
	double rx = nAtom.x()-cAtom.x();
	double ry = nAtom.y()-cAtom.y();
	double rz = nAtom.z()-cAtom.z();
	double norm = Helixer.norm(rx,ry,rz);
	rx /= norm;
	ry /= norm;
	rz /= norm;
 	// finding the rotation matrix
	double[][] rot = buildRotMatrix(rx,ry,rz,ang);
	double cmx=cAtom.x(),cmy=cAtom.y(),cmz=cAtom.z();
	
	for (int c=0 ; c<list.size() ; c++) {
		// Move atoms to C's origin
		double[] vec = {list.atomAt(c).x()-cmx,
		list.atomAt(c).y()-cmy,
		list.atomAt(c).z()-cmz};
		// rotating and translating 
		double[] tmp = Helixer.matTimesVec(rot,vec);
		list.atomAt(c).setXYZ(tmp[0]+cmx ,tmp[1]+cmy ,tmp[2]+cmz);
	}
}

public static double[][] buildRotMatrix(double x, double y, double z, double ang) {
	double cosang = Math.cos(ang);
	double sinang = Math.sin(ang);	
	double[][] rot = new double[3][3];
	rot[0][0] = cosang + (1-cosang)*x*x;
	rot[0][1] = (1-cosang)*x*y - sinang*z;
	rot[0][2] = (1-cosang)*x*z + sinang*y;
	rot[1][0] = (1-cosang)*x*y + sinang*z;	
	rot[1][1] = cosang + (1-cosang)*y*y;
	rot[1][2] = (1-cosang)*y*z - sinang*x;
	rot[2][0] = (1-cosang)*x*z - sinang*y;
	rot[2][1] = (1-cosang)*y*z + sinang*x;
	rot[2][2] = cosang + (1-cosang)*z*z;
	return rot;
}


public static Atom find(AtomList al , String name) {
	for (int c=0 ; c<al.size() ; c++)
		if (al.atomAt(c).name().equals(name))
             return al.atomAt(c);       
    return null;
}


}
