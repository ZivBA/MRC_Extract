package meshi.geometry;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.util.overlap.Overlap;

public class SmartOverlap {

public SmartOverlap() {}


/**
* Will move residues [r1,r2] in prot so that the residues in resList are adjacent.
* The structure of resList is:
* X1 Y1
* X2 Y2
* X3 Y3
*  ...
*
* {Y1,Y2...} must be in [r1,r2].
* The over lap will try to adjoin X1 and Y1, X2 and Y2 etc.
**/
public static void moveChunk(Protein protein , int r1 , int r2 , int[][] resList) { 
	AtomList ref = new AtomList();
	AtomList move = new AtomList();
	for (int c=0 ; c<resList.length ; c++) {
		Atom atom1 = protein.atoms().findAtomInList('G' , resList[c][0]);
		Atom atom2 = protein.atoms().findAtomInList('G' , resList[c][1]);
		if ((atom1!=null) && (atom2!=null)) {
			ref.add(protein.atoms().findAtomInList('G' , resList[c][0]));
			move.add(protein.atoms().findAtomInList('B' , resList[c][1]));
			ref.add(protein.atoms().findAtomInList('B' , resList[c][0]));
			move.add(protein.atoms().findAtomInList('G' , resList[c][1]));			
		}
		else {
			ref.add(protein.atoms().findAtomInList("CB" , resList[c][0]));
			move.add(protein.atoms().findAtomInList("CA" , resList[c][1]));
			ref.add(protein.atoms().findAtomInList("CA" , resList[c][0]));
			move.add(protein.atoms().findAtomInList("CB" , resList[c][1]));
		}
		 
	}
	makeOverlap(protein,r1,r2,ref,move);
}


public static void makeOverlap(Protein prot , int begin , int end , AtomList ref , AtomList al) {
	if (ref.size() != al.size())
		throw new RuntimeException("\n\nThe size of the reference list and the matching atom list must be the same. \n\n"); 
	
	// counting the number of atoms between start and end
	int total=0;
	for (int c=0 ; c<prot.atoms().size() ; c++) 
		if ((prot.atoms().atomAt(c).residueNumber()>=begin) &&
			(prot.atoms().atomAt(c).residueNumber()<=end)) 
			total++;
			
	// Building the coordinate arrays
	double[][] coRef = new double[3][total];
	double[][] coList = new double[3][total];
	int[] special = new int[ref.size()];
	int counter = 0;
	int specialCounter = 0;
	for (int c=0 ; c<prot.atoms().size() ; c++)
		if ((prot.atoms().atomAt(c).residueNumber()>=begin) &&
			(prot.atoms().atomAt(c).residueNumber()<=end)) {
				coList[0][counter] = prot.atoms().atomAt(c).x();
				coList[1][counter] = prot.atoms().atomAt(c).y();
				coList[2][counter] = prot.atoms().atomAt(c).z();
				int ind = al.findIndexInList(prot.atoms().atomAt(c));
				if (ind > -1) {
					coRef[0][counter] = ref.atomAt(ind).x();
					coRef[1][counter] = ref.atomAt(ind).y();
					coRef[2][counter] = ref.atomAt(ind).z();
					special[specialCounter] = counter;
					specialCounter++;
				}
				counter++;
			}
	 
	 // Finding the overlap
	 Overlap.rmsPartial(coRef, coList, special);
	 
	 // Puting the new coordinates into the atoms
	 counter = 0;
	 for (int c=0 ; c<prot.atoms().size() ; c++)
		if ((prot.atoms().atomAt(c).residueNumber()>=begin) &&
			(prot.atoms().atomAt(c).residueNumber()<=end)) {
				prot.atoms().atomAt(c).setXYZ(coList[0][counter],coList[1][counter],coList[2][counter]);
				counter++;
			}
}
	
}
