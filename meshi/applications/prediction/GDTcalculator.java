package meshi.applications.prediction;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.overlap.Overlap;

public class GDTcalculator   
{
	private static double[][] C1;
	private static double[][] C2;
	private static int[] convertingFromC1toResidueNum;
	private static int len;
	private static final int num_of_loops = 20;//the maximal number of iterations over the protein in finding the best superposition.


	public static double gdt(AtomList reference, AtomList protein1) {
		return gdt(reference, protein1, 1, 2, 4, 8);
	}

	public static double gdt(AtomList reference, AtomList protein1 , double cutoff1, double cutoff2, double cutoff3, double cutoff4) {
		double result = 0.0;
		int numBelow = 0;
		//creating two arrays (one for each given structure)
		//according to the atoms coordinates writen at the files.
		read_files(reference, protein1);

		// CUTOFF1        
		numBelow = findNumBelow(cutoff1);
		System.out.println(cutoff1 + " " + (numBelow*1.0/reference.CAFilter().size()));
		result += 0.25*numBelow*1.0/reference.CAFilter().size();

		// CUTOFF2        
		numBelow = findNumBelow(cutoff2,true);
		System.out.println(cutoff2 + " " + (numBelow*1.0/reference.CAFilter().size()));
		result += 0.25*numBelow*1.0/reference.CAFilter().size();

		// CUTOFF3        
		numBelow = findNumBelow(cutoff3);
		System.out.println(cutoff3 + " " + (numBelow*1.0/reference.CAFilter().size()));
		result += 0.25*numBelow*1.0/reference.CAFilter().size();

		// CUTOFF4        
		numBelow = findNumBelow(cutoff4);
		System.out.println(cutoff4 + " " + (numBelow*1.0/reference.CAFilter().size()));
		result += 0.25*numBelow*1.0/reference.CAFilter().size();

		return result;
	}//gdt

	public static int[] findAligningResidues(AtomList reference, AtomList protein1 , double cutoff) {
		//creating two arrays (one for each given structure)
		//according to the atoms coordinates writen at the files.
		read_files(reference, protein1);
		return findMatchingRes(cutoff);
	}

	public static void alignBySubset(AtomList reference, AtomList protein1 , double cutoff) {
		//creating two arrays (one for each given structure)
		//according to the atoms coordinates writen at the files.
		read_files(reference, protein1);
		int[] alignedRes = findAligningResidues(reference, protein1, cutoff);
		int[] alignIndices = new int[alignedRes.length];
		
		double[][] coordinatesRef = new double[3][protein1.size()];
		double[][] coordinatesMove = new double[3][protein1.size()];
		int alignedCounter=0;
		for (int c=0 ; c<protein1.size() ; c++) {
			coordinatesMove[0][c] = protein1.atomAt(c).x();
			coordinatesMove[1][c] = protein1.atomAt(c).y();
			coordinatesMove[2][c] = protein1.atomAt(c).z();

			Atom refAtom = reference.findAtomInList(protein1.atomAt(c).name(), protein1.atomAt(c).residueNumber());
			if (refAtom==null) {
//				System.out.println(" Warning (???): Could not find this atom in the reference:  " + protein1.atomAt(c));
				coordinatesRef[0][c] = 999;
				coordinatesRef[1][c] = 999;
				coordinatesRef[2][c] = 999;				
			}
			else {
				coordinatesRef[0][c] = refAtom.x();
				coordinatesRef[1][c] = refAtom.y();
				coordinatesRef[2][c] = refAtom.z();
				if (refAtom.name().equals("CA")) {
					boolean inAlignedList = false;
					for (int cc=0 ; cc<alignedRes.length ; cc++) {
						if (alignedRes[cc]==refAtom.residueNumber()) {
							inAlignedList=true;
						}
					}
					if (inAlignedList) {
						alignIndices[alignedCounter] = c;
						alignedCounter++;
					}				
				}
			}
		}
		
		System.out.println("The RMS on the subset is: " + Overlap.rmsPartialAltRMS(coordinatesRef, coordinatesMove, alignIndices));
		for (int c=0 ; c<protein1.size() ; c++) {
			protein1.atomAt(c).setXYZ(coordinatesMove[0][c], coordinatesMove[1][c], coordinatesMove[2][c]);
		}
	}
	
	//---------------------------------------------------------------------------------------
	public static int findNumBelow(double cutoff) {
		return findNumBelow(cutoff,false);
	}

	
	public static int findNumBelow(double cutoff, boolean toPrint) {
		int temp_best = 0;
		GDT_position bestPosition = null;
		for(int subInd=0; subInd<3; subInd++) {
			int sub = 3+2*subInd;  //the length of the initial sub-protein that we check.
			for(int pos=0; pos<(len+1-sub); pos++) {  //go over all initial fragments.
				GDT_position position = new GDT_position(sub, pos, cutoff, C1 , C2);
				position.find_best_conformation(num_of_loops);
				if (position.subset.length>temp_best) {
					temp_best = position.subset.length;
					bestPosition = position;
				}
			}
		}
		if (false && toPrint) {
			for (int c=0 ; c<bestPosition.subset.length ; c++)
				System.out.print(convertingFromC1toResidueNum[bestPosition.subset[c]]+ " ");
			System.out.println();
		}
		return temp_best;
	}

	private static int[] findMatchingRes(double cutoff) {
		int temp_best = 0;
		GDT_position bestPosition = null;
		for(int subInd=0; subInd<3; subInd++) {
			int sub = 3+2*subInd;  //the length of the initial sub-protein that we check.
			for(int pos=0; pos<(len+1-sub); pos++) {  //go over all initial fragments.
				GDT_position position = new GDT_position(sub, pos, cutoff, C1 , C2);
				position.find_best_conformation(num_of_loops);
				if (position.subset.length>temp_best) {
					temp_best = position.subset.length;
					bestPosition = position;
				}
			}
		}
		int[] result = new int[bestPosition.subset.length];
		for (int c=0 ; c<bestPosition.subset.length ; c++)
			result[c] = convertingFromC1toResidueNum[bestPosition.subset[c]];
		return result;
	}

	//-----------------------------------------------------------------------------------------------

	public static void read_files(AtomList reference, AtomList protein1) {
		AtomList al1,al2;
		int resnum,ind1,ind2;
		boolean found = false;

		len = 0;
		al1 = reference.CAFilter();
		al2 = protein1.CAFilter();
		for (ind2=0 ; ind2<al2.size() ; ind2++) {
			resnum = al2.atomAt(ind2).residueNumber();
			found = false;
			for (ind1=0 ; ((ind1<al1.size()) && !found)  ; ind1++) {
				if (resnum == al1.atomAt(ind1).residueNumber()) {
					if (!al2.atomAt(ind2).residueName().equals(al1.atomAt(ind1).residueName())) {
						System.out.print(al2.atomAt(ind2) + "\n" + al1.atomAt(ind1)+ "\n");
						throw new RuntimeException("The residue names in the above residues mismatch:\n\n");
					}
					len++;
					found = true;
				}
			}
		}
		
		C1 = new double[3][len];
		C2 = new double[3][len];
		convertingFromC1toResidueNum = new int[len];
		len = 0;
		for (ind2=0 ; ind2<al2.size() ; ind2++) {
			resnum = al2.atomAt(ind2).residueNumber();
			found = false;
			for (ind1=0 ; ((ind1<al1.size()) && !found)  ; ind1++) {
				if (resnum == al1.atomAt(ind1).residueNumber()) {
					convertingFromC1toResidueNum[len] = resnum;
					C1[0][len] =  al1.atomAt(ind1).x();
					C1[1][len] =  al1.atomAt(ind1).y();
					C1[2][len] =  al1.atomAt(ind1).z();
					C2[0][len] =  al2.atomAt(ind2).x();
					C2[1][len] =  al2.atomAt(ind2).y();
					C2[2][len] =  al2.atomAt(ind2).z();
					len++;
					found = true;
				}
			}
		}
		System.out.println("Found a match of " + len + " residues between reference and model, which are " + len*100.0/al1.size() + "% of the reference.");
	}//read_files

}//class 

// 

//# reference model file name is T0324_D2.pdb
//# initial model file name is T0324TS025_2-D2_res
//# Cutoff 1 is 1.0
//# Cutoff 2 is 2.0
//# Cutoff 3 is 4.0
//# Cutoff 4 is 8.0
//1.0 0.5538461538461539
//2.0 0.8
//4.0 0.9692307692307692
//8.0 1.0
//GDT: 0.8307692307692308      RMS: 1.7769353195720228




