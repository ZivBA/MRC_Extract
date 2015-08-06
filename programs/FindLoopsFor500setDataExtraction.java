package programs;

import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.molecularElements.Protein;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.dssp.CrystalContacts;
import meshi.util.dssp.DSSP;

/**
 * This application does not choose loops like 'FindLoopsFor500set' but just give the crystal contact
 * and exposure info. Later this data can be interpreted by MATLAB.
 * The output file has three columns:
 * Column 1: residue number.
 * Column 2: 0-if no crystal contact. 1-if crystal contact.
 * Column 3: Exposure [0.0 - 1.0]
 * 
 * run:
 * java FindLoopsFor500setDataExtraction <prot file name> <CryCo file name> <DSSP filename> <output file name> <crystal contact cutoff (4.0 is pretty good)>
 *
 */

public class FindLoopsFor500setDataExtraction implements Residues, AtomTypes { 

	public static void main(String[] args) {
		String protName = args[0].trim();
		String CryCoFileName = args[1].trim();
		String dsspFileName = args[2].trim();
		String outputFileName = args[3].trim();
		double crystalContactDis = (new Double(args[4].trim())).doubleValue();
		
		Protein prot =  new ExtendedAtomsProtein(protName,DO_NOT_ADD_ATOMS);
		CrystalContacts contacts = new CrystalContacts(CryCoFileName,crystalContactDis);
		DSSP dssp = new DSSP(dsspFileName);
		
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName));
			for (int resInd=0 ; resInd<prot.residues().size(); resInd++) {
				if (!prot.residueAt(resInd).dummy()) {
					int isCryCo = 0;
					if (contacts.isContact(prot.residueAt(resInd).number))
						isCryCo = 1;
					bw.write(prot.residueAt(resInd).number + " " + isCryCo + " " + 
							dssp.relACCofRes(prot.residueAt(resInd).number,' ') + "\n");
				}
			}
			bw.close();
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}	
	}
}
