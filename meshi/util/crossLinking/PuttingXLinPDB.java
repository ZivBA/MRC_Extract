package meshi.util.crossLinking;

import java.io.IOException;

import meshi.energy.AbstractEnergy;
import meshi.energy.TotalEnergy;
import meshi.energy.angle.AngleEnergy;
import meshi.energy.angle.AngleParametersList;
import meshi.energy.bond.BondEnergy;
import meshi.energy.bond.BondParametersList;
import meshi.energy.disConst.DisConstEnergy;
import meshi.energy.disConst.DisConstEnergyElement;
import meshi.energy.softExcludedVol.SoftExcludedVol;
import meshi.energy.softExcludedVol.SoftExcludedVolParametersList;
import meshi.energy.tether.TetherEnergy;
import meshi.geometry.AngleList;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.AtomPairList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LBFGS;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.file.MeshiWriter;

public class PuttingXLinPDB implements Residues, AtomTypes {
	
	private AtomList atomList = null;

	
	/**
	 * Constructor needs the atom list of the full complex.
	 * Note that the Atom instances in the list will be moved and frozen. 
	 */
	public PuttingXLinPDB(AtomList atomList) {
		this.atomList = atomList;
		atomList.moveCMtoOrigin();
//		try {
//			atomList.filter(new AtomList.NonHydrogen()).print(new MeshiWriter("fullComplex_centered.pdb"));
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		System.exit(1);
	}
	
	
	public void getElipsoide() {
		for (int c=0 ; c<atomList.size() ; c++) {
			if (atomList.atomAt(c).name().equals("CA")) {
				System.out.println(Math.pow(atomList.atomAt(c).x(),2) + " " + 
						Math.pow(atomList.atomAt(c).y(),2) + " " + 
						Math.pow(atomList.atomAt(c).z(),2));
			}
		}		
	}
	

	public int insideOutside(Crosslink xl) {
		Atom atom = atomList.findAtomInListReturningAtom("CB", xl.protName1(), xl.absPos1());
		int outsideScore = 0;
		if (atom!=null) {
			if (Math.pow(atom.x(),2) + 
					Math.pow(atom.y(),2) + 
					Math.pow(atom.z(),2) > Math.pow(68.0,2)) { // Outside
				outsideScore++;
			}
		}
		atom = atomList.findAtomInListReturningAtom("CB", xl.protName2(), xl.absPos2());
		if (atom!=null) {
			if (Math.pow(atom.x(),2) + 
					Math.pow(atom.y(),2) + 
					Math.pow(atom.z(),2) > Math.pow(68.0,2)) { // Outside
				outsideScore++;
			}
		}
		return outsideScore;
	}
	
	public void insideOutsideAllvec(CrosslinkVectorTRiC xlVec) {
		for (Crosslink xl : xlVec) {
			System.out.println(insideOutside(xl));
		}
	}
	
	
	public void refineXLatoms(Crosslink xl , int xlInd) {
		Atom atom1 = atomList.findAtomInListReturningAtom("CB", xl.protName1(), xl.absPos1());
		Atom atom2 = atomList.findAtomInListReturningAtom("CB", xl.protName2(), xl.absPos2());
		if ((atom1==null) | (atom2==null)) {
			System.out.println("SKIPPING: XL in unstructured part.");
			return;
		}
		Atom.resetNumberOfAtoms();
		AtomList al = new AtomList();
		for (int c=0 ; c<atomList.size() ; c++) {
			if ((atom1.distanceFrom(atomList.atomAt(c))<30.0) ||
					(atom2.distanceFrom(atomList.atomAt(c))<30.0)) {
				Atom newAtom = new Atom(atomList.atomAt(c).x(), 
						atomList.atomAt(c).y(), 
						atomList.atomAt(c).z(),
						atomList.atomAt(c).name(),
						atomList.atomAt(c).residueName(), 
						atomList.atomAt(c).residueNumber(),
						8);
				newAtom.setChain(atomList.atomAt(c).chain());
				al.add(newAtom);
			}
		}
		al.freeze();
		Protein twoLysProt = new Protein("aux2Lys.pdb",
				 new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		twoLysProt.defrost();
		stageXLatoms(xl,twoLysProt,atom1,atom2);
		al.add(twoLysProt.atoms());
		DistanceMatrix distanceMatrix = new DistanceMatrix(al, 5.5, 2.0, 4);  
	 	AtomPairList bondList = twoLysProt.bonds();
		bondList.renumber();
		BondEnergy bondEnergy = new BondEnergy(bondList, distanceMatrix, 
				new BondParametersList("C:/Users/Nir/Nir_Programs_eclipse/MESHI/meshicurrent/parameters/meshiPotential/bondEnergyParameters.dat") ,
				1.0);
		AngleList angleList = new AngleList(bondList, distanceMatrix);
		AngleEnergy angleEnergy = new AngleEnergy(angleList, distanceMatrix, 
				new AngleParametersList("C:/Users/Nir/Nir_Programs_eclipse/MESHI/meshicurrent/parameters/meshiPotential/angleEnergyParameters.dat") ,
				1.0);		
		SoftExcludedVol softEV = new SoftExcludedVol(distanceMatrix,
				new SoftExcludedVolParametersList("C:/Users/Nir/Nir_Programs_eclipse/MESHI/meshicurrent/parameters/meshiPotential/ExcludedVolumeParameters.dat"),
				1,
				1.0);
		DistanceList disList = new DistanceList();
		disList.add(twoLysProt.atoms().findAtomInList("CG", 1000).distance(twoLysProt.atoms().findAtomInList("CG", 1001)));
		disList.add(twoLysProt.atoms().findAtomInList("C", 1000).distance(
				al.findAtomInListReturningAtom("CA", xl.protName1(), xl.absPos1())));
		disList.add(twoLysProt.atoms().findAtomInList("C", 1001).distance(
				al.findAtomInListReturningAtom("CA", xl.protName2(), xl.absPos2())));
		DisConstEnergy disConstEnergy1 = new DisConstEnergy(disList,	1.0);	
		((DisConstEnergyElement) disConstEnergy1.elementsList().elementAt(0)).setTargetDis(1.6);
		((DisConstEnergyElement) disConstEnergy1.elementsList().elementAt(1)).setTargetDis(5.5);
		((DisConstEnergyElement) disConstEnergy1.elementsList().elementAt(2)).setTargetDis(5.5);
		DisConstEnergy disConstEnergy2 = new DisConstEnergy(disList,	1000.0);	
		((DisConstEnergyElement) disConstEnergy2.elementsList().elementAt(0)).setTargetDis(1.6);
		((DisConstEnergyElement) disConstEnergy2.elementsList().elementAt(1)).setTargetDis(5.5);
		((DisConstEnergyElement) disConstEnergy2.elementsList().elementAt(2)).setTargetDis(5.5);
		disConstEnergy2.off();
		AbstractEnergy[] energies = {bondEnergy , angleEnergy,  softEV, disConstEnergy1, disConstEnergy2}; 
		TotalEnergy energy = new TotalEnergy(al, distanceMatrix, energies);
		LBFGS minimizer = new LBFGS(energy, 0.0005, 10000, 100);
		try {
			System.out.println(minimizer.minimize());
			disConstEnergy2.on();
			System.out.println(minimizer.minimize());
		} catch (Exception e) {
			throw new RuntimeException("Minimization failed!!!");
		} 
		try {
			twoLysProt.atoms().print(new MeshiWriter("XL.pdb"));
		} catch (Exception e) {
			throw new RuntimeException("Writing failed!!!");
		} 
		
		for (int c=0 ; c<twoLysProt.atoms().size() ; c++) {
			Atom atom = twoLysProt.atoms().atomAt(c);
			Atom newAtom = new Atom(atom.x(), 
					atom.y(), 
					atom.z(),
					atom.name(),
					atom.residueName(), 
					atom.residueNumber() + 4*xlInd,
					8);
			System.out.println(newAtom);
		}
		
	}

	
	private void stageXLatoms(Crosslink xl, Protein twoLysProt, 
			Atom atom1, Atom atom2) {
		double x1 = atom1.x();
		double y1 = atom1.y();
		double z1 = atom1.z();		
		double x2 = atom2.x();
		double y2 = atom2.y();
		double z2 = atom2.z();
		int outside = insideOutside(xl);
		double rFactor = 0.7;
		if (outside>0) {
			rFactor = 1.3;
		}
		String[] names = {"O","C","CA","CB","CG","CG","CB","CA","C","O"};
		int[] resNums = {1000,1000,1000,1000,1000,1001,1001,1001,1001,1001};
		for (int frac=0 ; frac<10 ; frac++) {
			double x = rFactor*(x1 + (x2-x1)/9.0*frac + Math.random());
			double y = rFactor*(y1 + (y2-y1)/9.0*frac + Math.random());
			double z = rFactor*(z1 + (z2-z1)/9.0*frac + Math.random());
			Atom atom = twoLysProt.atoms().findAtomInList(names[frac], resNums[frac]);
			atom.setXYZ(x, y, z);
		}
	}
	

}
