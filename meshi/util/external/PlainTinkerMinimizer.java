package meshi.util.external;

import java.io.IOException;

import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.file.MeshiWriter;

public class PlainTinkerMinimizer implements Residues{
	
	
	public static void prepare(String doName , char unitName  , DunbrackLib lib) {
		TINKERandMESHIcomparison tink = new TINKERandMESHIcomparison(doName+".xyz", doName+".pdb");
		tink.fixBreaksInXYZ(doName+".tink.xyz");
		tink.setInactiveSuperPositionAtoms(unitName);
		tink.addStringToKeyFile(doName+".key", doName+".tink.key", tink.getInactiveString());
		tink.addStringToKeyFile(doName+".tink.key", doName+".tink.key", tink.getRestrainString());
		MESHIonTriC meshiOnTriC  = new MESHIonTriC(doName+".pdb", lib );	
		double[][] pp = meshiOnTriC.getNearestRotInfo();
		tink.addStringToKeyFile(doName+".tink.key", doName+".tink.key", tink.getTorsionRestrainString(pp));
	}
	
	
	
	public static void main(String[] args) {
		MeshiProgram.initRandom(0);
		if (args[0].equals("1")) {
			AtomList al = new AtomList(args[1]);
			al.setChain("A");
			try {
				al.print(new MeshiWriter("toMini.pdb"));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		else if (args[0].equals("2")) {
			String commandString = args[1];
			String doName = args[2];
			char unitName = args[3].trim().charAt(0);
			CommandList commands = new CommandList(commandString);
			DunbrackLib lib = new DunbrackLib(commands,0.99,100);
			PlainTinkerMinimizer.prepare(doName, unitName , lib);
		}
		else if (args[0].equals("3")) {
			Protein protOrig = new Protein(new AtomList(args[1]) , 
					new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
			Protein protMini = new Protein(new AtomList(args[2]) , 
					new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
			int treatingAtom = 0;
			AtomList newList = new AtomList();
			for (int resc=0 ; resc<protOrig.residues().size() ; resc++) {
				if (protOrig.residues().residueAt(resc).ca()!=null) {
					AtomList resList = protMini.atoms().getNextResidue(treatingAtom);
					treatingAtom += resList.size();
					for (int atc=0 ; atc<resList.size() ; atc++) {
						newList.add( new Atom(resList.atomAt(atc).x(),
								resList.atomAt(atc).y(), 
								resList.atomAt(atc).z(), 
								resList.atomAt(atc).name(), 
								resList.atomAt(atc).residueName(), 
								protOrig.residues().residueAt(resc).ca().residueNumber(),
								-1));
					}
				}
			}
			try {
				newList.print(new MeshiWriter("toMini.mini.pdb"));
			} catch (IOException e) {
				e.printStackTrace();
			}			
		}
	}

}
