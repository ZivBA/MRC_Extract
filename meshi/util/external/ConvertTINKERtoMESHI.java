package meshi.util.external;

import java.io.IOException;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.file.MeshiWriter;

public class ConvertTINKERtoMESHI {

	/**
	 * Run as: ConvertTINKERtoMESHI <TINKER pdb> <original PDB> <output PDB>
	 * 
	 * Will fix chain names and residue numbers of the TINKER output.
	 * 
	 * @param args
	 */
	
	public static void main(String[] args) {
		String tinkerPDB = args[0];
		String origPDB = args[1];
		String outputPDB = args[2];
		
		AtomList tinkerAL = new AtomList(tinkerPDB);
		AtomList origAL = new AtomList(origPDB);
		AtomList outputAL = new AtomList();
		
		int tinkerCounter = 0;
		int origCounter = 0;
		for ( ; tinkerCounter<tinkerAL.size() ; ) {
			AtomList tinkerRes = tinkerAL.getNextResidue(tinkerCounter);
			AtomList origRes = origAL.getNextResidue(origCounter);
			int tinkerResNum = tinkerRes.atomAt(0).residueNumber();			
			int origResNum = origRes.atomAt(0).residueNumber();			
			if (!tinkerRes.atomAt(0).residueName().equals(origRes.atomAt(0).residueName())) {
				throw new RuntimeException("Residue mismatch in tinker residue: " + tinkerResNum );
			}
			for (int c=0 ; c<tinkerRes.size() ; c++) {
				Atom newAtom = new Atom(tinkerRes.atomAt(c).x(),
						tinkerRes.atomAt(c).y(),
						tinkerRes.atomAt(c).z(),
						tinkerRes.atomAt(c).name(), tinkerRes.atomAt(c).residueName(), 
						origResNum, -1);
				newAtom.setChain(origRes.findAtomInList("CA", origResNum).chain());
				if (!newAtom.name().equals("OXT")) {
					outputAL.add(newAtom);
				}
			}
			tinkerCounter += tinkerRes.size();
			origCounter += origRes.size();
		}
		
		try {
			outputAL.print(new MeshiWriter(outputPDB));
		} catch (IOException e) {
			throw new RuntimeException("Could not write output file: " + outputPDB);
		}
	}
	
}
