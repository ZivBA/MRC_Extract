package meshi.applications.TriC;

import java.util.Vector;

public class VectorOfOrganisms {
	
	private Vector<ReadTriCofOrganism> vec;

	public VectorOfOrganisms() {
		vec = new Vector<ReadTriCofOrganism>();
	}
	
	public void add(ReadTriCofOrganism org) {
		if (!vec.isEmpty()) { 
			if (!checkIntegrity(org))
				throw new RuntimeException("Problem with the data integrity.");
			org.setStartInd(whereOrgHeaderStart(org));
		}
		else {
			org.setStartInd(0);
		}
		vec.add(org);
	}
	
	public String getProfile(char unit, int posInHeader) {
		String profile = "";
		for (ReadTriCofOrganism org : vec) {
			profile += org.getLetter(unit, posInHeader);
		}
		return profile;
	}
	
	private boolean checkIntegrity(ReadTriCofOrganism org) {
		if (whereOrgHeaderStart(org)<0)
			return false;
		else
			return true;
	}
	
	private int whereOrgHeaderStart(ReadTriCofOrganism org) {
		return vec.get(0).headerSeq.replaceAll("-", "").indexOf(org.headerSeq.replaceAll("-", ""));
	}
	

}
