package meshi.util.TRIC.EM;

import java.io.IOException;

import meshi.applications.TriC.TricAlignment;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.parameters.Residues;
import meshi.util.TRIC.PutUnitInAnyTopPosition;
import meshi.util.file.MeshiWriter;

public class DockThermosomIntoEM_V2 implements Residues {

	public static boolean[] makeToTake() {
		boolean[] toTake = new boolean[522];
		for (int c=0 ; c<toTake.length ; c++) {
			toTake[c]=false;
		}
		// Doing EQUATORIAL domain
		for (int c=21 ; c<=37 ; c++) {
			toTake[c]=true;
		}
		for (int c=95 ; c<=116 ; c++) {
			toTake[c]=true;
		}
		for (int c=120 ; c<=140 ; c++) {
			toTake[c]=true;
		}
		for (int c=410 ; c<=424 ; c++) {
			toTake[c]=true;
		}
		for (int c=429 ; c<=452 ; c++) {
			toTake[c]=true;
		}
		for (int c=456 ; c<=468 ; c++) {
			toTake[c]=true;
		}
		for (int c=496 ; c<=513 ; c++) {
			toTake[c]=true;
		}
		// Doing MIDDLE domain
		for (int c=147 ; c<=159 ; c++) {
			toTake[c]=true;
		}
		for (int c=167 ; c<=182 ; c++) {
			toTake[c]=true;
		}
		for (int c=197 ; c<=202 ; c++) {
			toTake[c]=true;
		}
		for (int c=210 ; c<=213 ; c++) {
			toTake[c]=true;
		}
		for (int c=368 ; c<=374 ; c++) {
			toTake[c]=true;
		}		
		for (int c=380 ; c<=402 ; c++) {
			toTake[c]=true;
		}
		// Doing APICAL domain
		for (int c=216 ; c<=219 ; c++) {
			toTake[c]=true;
		}
		for (int c=230 ; c<=239 ; c++) {
			toTake[c]=true;
		}
		for (int c=277 ; c<=282 ; c++) {
			toTake[c]=true;
		}
		for (int c=287 ; c<=290 ; c++) {
			toTake[c]=true;
		}
		for (int c=296 ; c<=304 ; c++) {
			toTake[c]=true;
		}
		for (int c=308 ; c<=311 ; c++) {
			toTake[c]=true;
		}
		for (int c=315 ; c<=324 ; c++) {
			toTake[c]=true;
		}
		for (int c=341 ; c<=349 ; c++) {
			toTake[c]=true;
		}
		for (int c=357 ; c<=364 ; c++) {
			toTake[c]=true;
		}
		// Doing the LID
		for (int c=252 ; c<=256 ; c++) {
			toTake[c]=true;
		}		
		for (int c=260 ; c<=272 ; c++) {
			toTake[c]=true;
		}
		return toTake;
	}

	

	public static void main(String[] args) {

		boolean[] toTake = makeToTake();
		int[][][] domainParsing = {{{17 , 144}, {403 , 520}},
				{{147 , 214} , {366 , 402}},
				{{215 , 244} , {275 , 365}},
				{{250 , 273}}};
		TricAlignment tricAlignment = new TricAlignment();
		
		// Looping on positions.
		for (int posC = 0 ; posC<8 ; posC++) {
			AtomList[] refs = new AtomList[4];
			refs[0] = new AtomList("../Maps/SetupE_pos_"+posC+"_full.pdb"); // Equatorial
			refs[1] = new AtomList("../Maps/SetupC_pos_"+posC+"_full.pdb"); // Middle
			refs[2] = new AtomList("../Maps/SetupA_pos_"+posC+"_full.pdb"); // Apical
			refs[3] = new AtomList("../Maps/SetupG_pos_"+posC+"_full.pdb"); // Lid
			// Looping on units
	        String units = "ABGDEHQZ";
			for (int unitC = 0 ; unitC<8 ; unitC++) {
				Atom.resetNumberOfAtoms();
				AtomList unit = new AtomList("../Maps/Atemp_HM_"+units.charAt(unitC)+".scwrl.pdb");
				AtomList writeUnitRaw = new AtomList();
				// Adding the atoms
				for (int domainC=0 ; domainC<domainParsing.length ; domainC++) {
					PutUnitInAnyTopPosition.alignByEquatorialDomains(refs[domainC], 'I', unit, units.charAt(unitC));
					for (int domainPart=0 ; domainPart<domainParsing[domainC].length ; domainPart++) {
						for (int atomC=0 ; atomC<unit.size() ; atomC++) {
							int resNumInThermo = tricAlignment.getNewResNum(units.charAt(unitC), unit.atomAt(atomC).residueNumber(), 'I');
							if ((resNumInThermo>=domainParsing[domainC][domainPart][0]) &&
									(resNumInThermo<=domainParsing[domainC][domainPart][1])) {
								if (toTake[resNumInThermo]) {
									Atom atom = new Atom(unit.atomAt(atomC).x(),
											unit.atomAt(atomC).y(),
											unit.atomAt(atomC).z(),
											unit.atomAt(atomC).name(),
											unit.atomAt(atomC).residueName(),
											unit.atomAt(atomC).residueNumber(),
											-1);
									atom.resetEnergy();
									atom.addEnergy(resNumInThermo);
									writeUnitRaw.add(atom);
								}
							}							
						}
					}
				}
				// Sorting the atom list
				AtomList writeUnit = new AtomList();
				for (int resC=0 ; resC<600 ; resC++) {
					for (int atomC=0 ; atomC<writeUnitRaw.size() ; atomC++) {
						if (writeUnitRaw.atomAt(atomC).residueNumber()==resC) {
							writeUnit.add(writeUnitRaw.atomAt(atomC));
						}
					}					
				}
				// Writing to disk
				try {
					writeUnit.print(new MeshiWriter("../Nir_subunits/AA_SS_"+posC+"_"+units.charAt(unitC)+".pdb"));
				} catch (IOException e) {
					throw new RuntimeException("Could not write file.");
				}	
			}
		}
	} // Of main
	
	
}
