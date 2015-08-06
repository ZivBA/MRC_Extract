package meshi.energy.torsionSpaceMinimization;

import java.util.Vector;

import meshi.geometry.Angle;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionList;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;

public class TorsionMinimizationTree {
	
	private Protein prot;
	private DistanceMatrix dm;
	private Vector<TorsionMinimizationElement> tree;
	private TempCalcForthAtom tempCalcForthAtom = new TempCalcForthAtom(); 
	
	public TorsionMinimizationTree(Protein prot, DistanceMatrix dm){
		this.prot = prot;
		this.dm = dm;
		for (int c=0 ; c<prot.atoms().size(); c++)
			prot.atoms().atomAt(c).wasVisited = false;
		tree = new Vector<TorsionMinimizationElement>();
		buildTree();
	}
	
	
	private void buildTree() {
		for (int c=0 ; c<prot.residues().size() ; c++) 
			if (!prot.residues().residueAt(c).dummy()) 
				if ((c==0) ||  prot.residues().residueAt(c-1).dummy() ||
						((prot.residues().residueAt(c-1).number+1) != prot.residues().residueAt(c).number)) {
					int firstResNum = prot.residues().residueAt(c).number;
/*					Atom dummyFirstC = getDummyFirstC(firstResNum, 
							prot.atoms().findAtomInList("C",firstResNum),
							prot.atoms().findAtomInList("CA",firstResNum),
							prot.atoms().findAtomInList("N",firstResNum));
					prot.atoms().findAtomInList("CA",firstResNum).wasVisited = true;
					prot.atoms().findAtomInList("N",firstResNum).wasVisited = true;
					dummyFirstC.wasVisited = true;
					DFS(dummyFirstC,prot.atoms().findAtomInList("N",firstResNum),
							prot.atoms().findAtomInList("CA",firstResNum));              */
					prot.atoms().findAtomInList("C",firstResNum).wasVisited = true;
					prot.atoms().findAtomInList("CA",firstResNum).wasVisited = true;
					prot.atoms().findAtomInList("N",firstResNum).wasVisited = true;
					DFS(prot.atoms().findAtomInList("N",firstResNum),
							prot.atoms().findAtomInList("CA",firstResNum),
							prot.atoms().findAtomInList("C",firstResNum));					
				}
	}
	

	public void buildProtein() {
		for (TorsionMinimizationElement tme : tree) {
			tme.buildTorsionInProt(tempCalcForthAtom);
		}
	}
	
	
	private TorsionMinimizationElement DFS(Atom a1, Atom a2, Atom a3) {
		int whereMainTorsion = whereMainTorsion(a3.bonded(), a1, a2, a3);
		if (whereMainTorsion==-1) { // Leaf...
			return null;
		}
		a3.bonded().atomAt(whereMainTorsion).wasVisited = true;
		TorsionMinimizationElement tme = new TorsionMinimizationElement(a1,a2,a3,a3.bonded().atomAt(whereMainTorsion),dm);
		tree.add(tme);		
		tme.setSonMain(DFS(a2,a3,a3.bonded().atomAt(whereMainTorsion)));
		for (int c=0 ; c<a3.bonded().size(); c++) {  // Other branches... 
			if (!a3.bonded().atomAt(c).wasVisited) {
				a3.bonded().atomAt(c).wasVisited = true;
				tme.addBranch(a3.bonded().atomAt(c), DFS(a2,a3,a3.bonded().atomAt(c)));
			}
		}
		return tme;
	}
	
	
	private int whereMainTorsion(AtomList list, Atom a1, Atom a2, Atom a3) {
		/*
		// If ever I could use this first dummy atom...
		// Special cases:
		if (a1.type==-1) { // We are in the first "N" and "CA" of the chain
			for (int c=0; c<list.size() ; c++) 
				if (!list.atomAt(c).wasVisited)
					if (list.atomAt(c).name().equals("C"))
						return c;
			return -1;
		}
		*/
		for (int c=0; c<list.size() ; c++) 
			if (!list.atomAt(c).wasVisited) {
				Torsion tor = new Torsion(new Angle(a1, a2, a3, dm, false), 
						new Angle(a2, a3, list.atomAt(c), dm, false), dm);
				if (TorsionList.isNameInIUPAC(tor)) 
					return c;
			}
		for (int c=0; c<list.size() ; c++) 
			if (!list.atomAt(c).wasVisited)
				return c;		
		return -1;
	}

/*
	// If ever I could use this first dummy atom...
	private Atom getDummyFirstC(int resNum, Atom C, Atom CA, Atom N) {
		double[] coors = tempCalcForthAtom.getCoors(1.33, (Math.PI/180)*116.7 , (Math.PI/180)*(-60), 
				C, CA, N);
		return new Atom(coors[0], coors[1], coors[2],
				   "C", "UNK", resNum-1, -1);
	}
*/
	
	/**
	 * updating the derivatives in the tree.
	 **/
	public void calcDerivatives() {
		for (int c=tree.size()-1 ; c>=0 ; c--){
			tree.get(c).calcDericatives();
		}
	}
	
	
	// Return all the elements around torsions with IUPAC naming except for elements around OMEGA torsions.
	public double[][] getCoors() {
		int numberOfcoors = 0;
		for (int c=0 ; c<tree.size() ; c++) 
			if (tree.get(c).isElementAroundSignificantTorsion())
				numberOfcoors++;
		double[][] coors = new double[numberOfcoors][];
		numberOfcoors = 0;
		for (int c=0 ; c<tree.size() ; c++) {
			if (tree.get(c).isElementAroundSignificantTorsion()) {
				coors[numberOfcoors] = tree.get(c).getCoors();
				numberOfcoors++;
			}
		}
		return coors;
	}
	
		
}
