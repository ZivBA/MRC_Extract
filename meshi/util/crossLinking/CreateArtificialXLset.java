package meshi.util.crossLinking;

import meshi.molecularElements.AtomList;

public class CreateArtificialXLset {

	private AtomList atomList = null;

	/**
	 * The purpose of this class is to create a cross-link set of any size from a certain complex.
	 */
	public CreateArtificialXLset(AtomList atomList) {
		this.atomList = atomList.filter(new AtomList.KCA_Filter());
	}
	
	public String createSet(int setSize , double maxXLdistance) {
		String outString = "";
		int maxNumberOfTries = 10;
		for (int soFarFound = 0; soFarFound<setSize ; ) {
			int atom1 = (int) (Math.random()*atomList.size());
			for (int tries = 0; tries<maxNumberOfTries ; tries++) {
				int atom2 = (int) (Math.random()*atomList.size());
				if (atomList.atomAt(atom1).distanceFrom(atomList.atomAt(atom2))<maxXLdistance) {
					if (!(atomList.atomAt(atom1).chain().equals(atomList.atomAt(atom2).chain()) &&
							atomList.atomAt(atom1).residueNumber() == atomList.atomAt(atom2).residueNumber())) {
						Crosslink xl = new Crosslink(-1, "K", "K", 0, 0, atomList.atomAt(atom1).chain(), atomList.atomAt(atom2).chain(), "N/A", atomList.atomAt(atom1).residueNumber(), atomList.atomAt(atom2).residueNumber());
						outString += (xl.toStringKalismanFormat() + "\n");
						soFarFound++;
						break;
					}
				}
			}
		}
		return outString;
	}
	

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		CreateArtificialXLset setMaker = new CreateArtificialXLset(new AtomList("One_Ring_Model_OMS.pdb"));
		System.out.println(setMaker.createSet(1000, 28.0));
	}

}
