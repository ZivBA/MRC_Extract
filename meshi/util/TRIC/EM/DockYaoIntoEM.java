package meshi.util.TRIC.EM;

import java.io.IOException;

import meshi.applications.TriC.TricAlignment;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.TRIC.PutUnitInAnyTopPosition;
import meshi.util.file.MeshiWriter;

/**
 * This scripts write the secondary structured parts of Michael's models on Yao's PDB.
 * 
 * @author Nir
 *
 */

public class DockYaoIntoEM {
	
	public static void main(String[] args) {

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
		for (int c=380 ; c<=402 ; c++) {
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
		
		// Putting the different units in different positions
		// Looping on positions.
		TricAlignment tricAlignment = new TricAlignment();
		for (int yaoC = 0 ; yaoC<8 ; yaoC++) {
			Atom.resetNumberOfAtoms();
	        String posString;
	        switch (yaoC) {
	        	case 0:  posString = "A_pos_00";       break;
	            case 1:  posString = "Z_pos_01";       break;
	            case 2:  posString = "B_pos_02";      break;
	            case 3:  posString = "G_pos_03";         break;
	            case 4:  posString = "Q_pos_04";         break;
	            case 5:  posString = "D_pos_05";           break;
	            case 6:  posString = "E_pos_06";          break;
	            case 7:  posString = "H_pos_07";          break;
	            default: throw new RuntimeException("Invalid position");
	        }	
			// Looping on units
	        String units = "ABGDEHQZ";
			for (int unitC = 0 ; unitC<8 ; unitC++) {
				Atom.resetNumberOfAtoms();
				AtomList anyUnitFull = (new AtomList(units.charAt(unitC)+"=3IYG_"+posString+".scwrl4.pdb")).noOXTFilter().filter(new AtomList.NonHydrogen());

				// Trimming to the secondary structure
				AtomList anyUnit = new AtomList();
				for (int c=0 ; c<anyUnitFull.size() ; c++) {
					int resNum = anyUnitFull.atomAt(c).residueNumber();
					int resNumInThermo = tricAlignment.getNewResNum(units.charAt(unitC), resNum, 'I');
					if (toTake[resNumInThermo]) {
						anyUnitFull.atomAt(c).resetEnergy();
						anyUnitFull.atomAt(c).addEnergy(resNumInThermo);
						anyUnit.add(anyUnitFull.atomAt(c));
					}
				}				
				anyUnit.setChain(units.charAt(unitC)+"");
				
				// Writing to disk
				try {
					anyUnit.print(new MeshiWriter("AA_SS_"+yaoC+"_"+units.charAt(unitC)+".pdb"));
				} catch (IOException e) {
					throw new RuntimeException("Could not write file.");
				}				
			}
		}
		
	} // Of main
	
	
	
	
}
