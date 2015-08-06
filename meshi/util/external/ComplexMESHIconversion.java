package meshi.util.external;

import java.io.IOException;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.Residues;
import meshi.util.file.MeshiWriter;

public class ComplexMESHIconversion implements Residues {

	
	public static Protein complex2meshi(String complexPDB) {
		AtomList complex = new AtomList(complexPDB);
		complex.defrost();
		return complex2meshi(complex);
	}
	
	public static Protein complex2meshi(AtomList complex) {
		AtomList meshiList = new AtomList();
		for (int c=0 ; c<complex.size() ; c++) {
			Atom meshiAtom = new Atom(complex.atomAt(c).x(),
					complex.atomAt(c).y(),
					complex.atomAt(c).z(),
					complex.atomAt(c).name(),
					complex.atomAt(c).residueName(),
					complex2meshi_resnum(complex.atomAt(c).chain(), complex.atomAt(c).residueNumber()),
					-1);
			meshiAtom.setChain("");
			if (complex.atomAt(c).frozen()) {
				meshiAtom.freeze();
			}
			meshiList.add( meshiAtom );
		}
		return new Protein(meshiList , new ResidueExtendedAtoms(ADD_ATOMS));
	}

	
	public static AtomList MEHSI2complex(Protein meshiProt) {
		AtomList complex = new AtomList();
		for (int c=0 ; c<meshiProt.atoms().size() ; c++) {
			Atom complexAtom = new Atom(meshiProt.atoms().atomAt(c).x(),
					meshiProt.atoms().atomAt(c).y(),
					meshiProt.atoms().atomAt(c).z(),
					meshiProt.atoms().atomAt(c).name(),
					meshiProt.atoms().atomAt(c).residueName(),
					meshi2complex_resnum(meshiProt.atoms().atomAt(c).residueNumber()),
					-1);
			complexAtom.setChain(meshi2complex_chain(meshiProt.atoms().atomAt(c).residueNumber()));
			complex.add(complexAtom);
		}
		return complex;
	}

	public static void writeMEHSI2complex(Protein meshiProt, String complexPDB) {
		AtomList complex = MEHSI2complex(meshiProt);
		try {
			complex.filter(new AtomList.NonHydrogen()).print(new MeshiWriter(complexPDB));
		} catch (IOException e) {
			throw new RuntimeException("Could not write file: " + complexPDB);
		}
	}

	
	public static int findAtomInOrigAtomList(AtomList alOrig, Atom atomInMESHI) {
		return alOrig.findAtomInList(atomInMESHI.name(),
				meshi2complex_chain(atomInMESHI.residueNumber()),
				meshi2complex_resnum(atomInMESHI.residueNumber()));
	}
	
	
	 // All the conversion smartness is in these three methods.
	/*
	 * Returns the MESHI residue number
	 */
	private static int complex2meshi_resnum(String chain, int resNum) {
		if (chain.equals("A")) {
			return resNum;
		}
		if (chain.equals("B")) {
			return 1000+resNum;
		}
		throw new RuntimeException("Currently can only handle a complex with chain A and B");
	}

	/*
	 * Returns the complex chain
	 */
	private static String meshi2complex_chain(int resNum) {
		if (resNum<1000) {
			return "A"; 
		}
		else {
			return "B";
		}
	}

	/*
	 * Returns the complex resnum
	 */
	private static int meshi2complex_resnum(int resNum) {
		if (resNum<1000) {
			return resNum; 
		}
		else {
			return resNum-1000;
		}
	}
	
}
