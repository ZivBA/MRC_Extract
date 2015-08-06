package programs;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.molecularElements.Protein;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.dssp.CrystalContacts;
import meshi.util.dssp.DSSP;

public class FindLoopsFor500Set implements Residues, AtomTypes {

	public static final int loopL = 8; 
	public static final int loopOverlap = 3; 
	public static final double crystalContactDis = 4.0; 
	public static final String protPath = "proteins/";
	public static final int minNumberOfCoils = 3;
	public static final int minNumberOfExposed = 2;	 
	public static final double exposureTH = 0.6;	 

	public static void main(String[] args) {
		String protName = args[0].trim();
		String CryCoFileName = args[1].trim();
		String dsspFileName = args[2].trim();
		
		Protein prot =  new ExtendedAtomsProtein(protPath+protName,DO_NOT_ADD_ATOMS);
		CrystalContacts contacts = new CrystalContacts(CryCoFileName,crystalContactDis);
		DSSP dssp = new DSSP(dsspFileName);

//		for (int c=0 ; c<contacts.getContacts().length ; c++)
	//		System.out.println(contacts.getContacts()[c]);
		for (int resInd=0 ; resInd<prot.residues().size(); resInd++) {
			int resNum = prot.residues().residueAt(resInd).number;
			if (loopResiduesExists(prot,resNum))
				if (loopResiduesNonContact(contacts,resNum))
					if (enoughCoil(dssp,resNum))
						if (enoughExposure(dssp,resNum))
							System.out.println(protName +  " " + resNum + " " + (resNum+loopL-1));
		}
	}

	// The loop and overlapping residues must exist
	private static boolean loopResiduesExists(Protein prot, int resNum) {
		for (int c=resNum-loopOverlap ; c<resNum+loopL+loopOverlap ; c++)
			if ((prot.residue(c)==null) || (prot.residue(c).dummy()))
				return false;
		return true;
	}
	

	// The loop and overlapping residues must not be in a crystal contact
	private static  boolean loopResiduesNonContact(CrystalContacts contacts, int resNum) {
		for (int c=resNum ; c<resNum+loopL ; c++)
		//for (int c=resNum-loopOverlap ; c<resNum+loopL+loopOverlap ; c++)
			if (contacts.isContact(c))
				return false;
		return true;		
	}
	
	// The loop and overlapping residues must not be in a crystal contact
	private static boolean enoughCoil(DSSP dssp, int resNum){
		int coilCounter = 0;
		for (int c=resNum ; c<resNum+loopL ; c++)
			if (((dssp.SSofRes(c,' ')!='H') && (dssp.SSofRes(c,' ')!='E')) || 
					(dssp.SSofRes(c,' ')!=dssp.SSofRes(c-1,' ')))	
				coilCounter++;
		if (coilCounter<minNumberOfCoils)
			return false;
		return true;				
	}
	
	// The loop and overlapping residues must not be in a crystal contact
	private static boolean enoughExposure(DSSP dssp, int resNum) {
		int exposedCounter = 0;
		for (int c=resNum ; c<resNum+loopL ; c++)
			if (dssp.relACCofRes(c,' ')>exposureTH)	
				exposedCounter++;
		if (exposedCounter<minNumberOfExposed)
			return false;
		return true;						
	}
	
}
