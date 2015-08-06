package meshi.util.external;

import java.io.IOException;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.file.MeshiWriter;
import meshi.util.overlap.Overlap;

public class ConvertNCNtoMESHI {

	/**
	 * Run as: ConvertNCNtoMESHI <NCN pdb> <original PDB> <output PDB>
	 * 
	 * Will translate-rotate each residue so that the backbone N,CA,C match is optimal.
	 * 
	 * @param args
	 */
	
	public static void main(String[] args) {
		String ncnPDB = args[0];
		String origPDB = args[1];
		String outputPDB = args[2];
		
		AtomList ncnAL = new AtomList(ncnPDB);
		AtomList origAL = new AtomList(origPDB);
		AtomList outputAL = new AtomList();
		
		int ncnCounter = 0;
		int origCounter = 0;
		for ( ; ncnCounter<ncnAL.size() ; ) {
			AtomList ncnRes = ncnAL.getNextResidue(ncnCounter);
			AtomList origRes = origAL.getNextResidue(origCounter);
			int ncnResNum = ncnRes.atomAt(0).residueNumber();			
			int origResNum = origRes.atomAt(0).residueNumber();			
			if (!ncnRes.atomAt(0).residueName().equals(origRes.atomAt(0).residueName())) {
				throw new RuntimeException("Residue mismatch in NCN residue: " + ncnResNum );
			}
			// Making the superpos arrays
			double[][] origAr = new double[3][3];
			double[][] ncnAr = new double[3][3+ncnRes.size()];
			Atom origN = origRes.findAtomInList("N", origResNum);
			origAr[0][0] = origN.x();
			origAr[1][0] = origN.y();
			origAr[2][0] = origN.z();
			Atom origCA = origRes.findAtomInList("CA", origResNum);
			origAr[0][1] = origCA.x();
			origAr[1][1] = origCA.y();
			origAr[2][1] = origCA.z();
			Atom origC = origRes.findAtomInList("C", origResNum);
			origAr[0][2] = origC.x();
			origAr[1][2] = origC.y();
			origAr[2][2] = origC.z();
			Atom ncnN = ncnRes.findAtomInList("N", ncnResNum);
			ncnAr[0][0] = ncnN.x();
			ncnAr[1][0] = ncnN.y();
			ncnAr[2][0] = ncnN.z();
			Atom ncnCA = ncnRes.findAtomInList("CA", ncnResNum);
			ncnAr[0][1] = ncnCA.x();
			ncnAr[1][1] = ncnCA.y();
			ncnAr[2][1] = ncnCA.z();
			Atom ncnC = ncnRes.findAtomInList("C", ncnResNum);
			ncnAr[0][2] = ncnC.x();
			ncnAr[1][2] = ncnC.y();
			ncnAr[2][2] = ncnC.z();
			for (int c=0 ; c<ncnRes.size() ; c++) {
				Atom ncnAtom = ncnRes.atomAt(c);
				ncnAr[0][3+c] = ncnAtom.x();
				ncnAr[1][3+c] = ncnAtom.y();
				ncnAr[2][3+c] = ncnAtom.z();				
			}
			int[] tmp = {0,1,2};
			Overlap.rmsPartialAltRMS(origAr, ncnAr, tmp);
			for (int c=0 ; c<ncnRes.size() ; c++) {
				Atom newAtom = new Atom(ncnAr[0][3+c], ncnAr[1][3+c], ncnAr[2][3+c],
						ncnRes.atomAt(c).name(), ncnRes.atomAt(c).residueName(), 
						origResNum, -1);
				newAtom.setChain(origRes.findAtomInList("CA", origResNum).chain());
				outputAL.add(newAtom);
			}
			ncnCounter += ncnRes.size();
			origCounter += origRes.size();
		}
		
		try {
			outputAL.print(new MeshiWriter(outputPDB));
		} catch (IOException e) {
			throw new RuntimeException("Could not write output file: " + outputPDB);
		}
	}
}
