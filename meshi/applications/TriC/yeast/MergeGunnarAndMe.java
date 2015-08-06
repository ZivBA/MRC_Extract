package meshi.applications.TriC.yeast;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.file.File2StringArray;

public class MergeGunnarAndMe {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String nir =    "C:\\Users\\Nir\\TRiC\\Crystallography\\ADP\\ADP_positions\\Put_K\\refine_17_fixed_ADP_withK.pdb";
		String gunnar = "C:\\Users\\Nir\\TRiC\\Crystallography\\ADP\\ADP_positions\\Put_K\\refine_17_modified.pdb";
		String out =    "C:\\Users\\Nir\\TRiC\\Crystallography\\ADP\\ADP_positions\\Put_K\\refine_17_fixed_ADP_withK_with_SF.pdb";
		AtomList fullList = (new AtomList(nir)).filter(new AtomList.NonHydrogen()).noOXTFilter();
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(out));
			String[] gFile = File2StringArray.f2a(gunnar);
			for (int c=0 ; c<gFile.length ; c++) {
				if ((gFile[c].length()<5) || !gFile[c].substring(0, 5).equals("ATOM ")) {
					bw.write(gFile[c] + "\n");
				}
				else if (gFile[c].contains("OXT")) {
					bw.write(gFile[c] + "\n");
				}
				else {
					Atom atom = fullList.findAtomInListReturningAtom(gFile[c].substring(12,16).trim(), gFile[c].substring(21,22), Integer.parseInt(gFile[c].substring(23,26).trim()));
					String atomS = atom.toString();
					if (gFile[c].substring(13, 16).equals("C5*") |
							gFile[c].substring(13, 16).equals("C4*") |
							gFile[c].substring(13, 16).equals("O4*") |
							gFile[c].substring(13, 16).equals("C3*") |
							gFile[c].substring(13, 16).equals("O3*") |
							gFile[c].substring(13, 16).equals("C2*") |
							gFile[c].substring(13, 16).equals("O2*") |
							gFile[c].substring(13, 16).equals("C1*") |
							gFile[c].substring(13, 16).equals("N9 ") |
							gFile[c].substring(13, 16).equals("C8 ") |
							gFile[c].substring(13, 16).equals("N7 ") |
							gFile[c].substring(13, 16).equals("C5 ") |
							gFile[c].substring(13, 16).equals("C6 ") |
							gFile[c].substring(13, 16).equals("N6 ") |
							gFile[c].substring(13, 16).equals("N1 ") |
							gFile[c].substring(13, 16).equals("C2 ") |
							gFile[c].substring(13, 16).equals("N3 ") |
							gFile[c].substring(13, 16).equals("C4 ") ) {
						bw.write(gFile[c] + "\n");
						System.out.print("kept from original: " + gFile[c] + "\n");
					}
					else {
						bw.write(gFile[c].substring(0, 26) + atomS.substring(26,54) + gFile[c].substring(54) + "\n");
						System.out.print(gFile[c].substring(0, 26) + atomS.substring(26,54) + gFile[c].substring(54) + "\n");
					}
				}
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException(e);
		}
		


	}

}
