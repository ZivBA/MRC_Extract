package meshi.util.crossLinking;

import meshi.molecularElements.AtomList;
import meshi.util.TRIC.PutUnitInAnyTopPositionYeast;

public class CrosslinkVectorTinghe extends CrosslinkVector {

	private static final long serialVersionUID = 1L;

	public CrosslinkVectorTinghe() {
		// TODO Auto-generated constructor stub
	}

	public CrosslinkVectorTinghe(String filename , int fileType) {
		super(filename ,  fileType);
		// TODO Auto-generated constructor stub
	}

	/**
	 * In Tinghe the stochiometry of all the subunits is 1.
	 */
	public String identicalChains(String protName) {
		return protName;
	}

	public CrosslinkVector createEmptyVector() {
		return new CrosslinkVectorTinghe();
	}	
	
	public static void main(String[] args) {
		
		CrosslinkVectorTinghe xlVec = new CrosslinkVectorTinghe("Excel_20120815_TWu_Pugilisi_HCD_Rigi_Red_02-high-confidence.txt",1);
		AtomList fullComplex = new AtomList("Model_ABE.pdb");
		for (Crosslink xl : xlVec) {
			String str = "N/A";
			if ((fullComplex.findAtomInListReturningAtom("CA", xl.protName1(), xl.absPos1()) != null) &
					(fullComplex.findAtomInListReturningAtom("CA", xl.protName2(), xl.absPos2()) != null)) {
				double disAA = fullComplex.findAtomInListReturningAtom("CA", "A", xl.absPos1()).distanceFrom(
						fullComplex.findAtomInListReturningAtom("CA", "A", xl.absPos2()) );
				double disAB = fullComplex.findAtomInListReturningAtom("CA", "A", xl.absPos1()).distanceFrom(
						fullComplex.findAtomInListReturningAtom("CA", "B", xl.absPos2()) );
				double disBA = fullComplex.findAtomInListReturningAtom("CA", "A", xl.absPos1()).distanceFrom(
						fullComplex.findAtomInListReturningAtom("CA", "E", xl.absPos2()) );
				double disAE = fullComplex.findAtomInListReturningAtom("CA", "E", xl.absPos1()).distanceFrom(
						fullComplex.findAtomInListReturningAtom("CA", "A", xl.absPos2()) );
				str = (((int) Math.min(Math.min(Math.min(disAA,disAB),disBA),disAE)*10)/10.0) + " " + (((int) disAA*10)/10.0) + " " + (((int) disAE*10)/10.0) + " " + (((int) disAB*10)/10.0) + " " + (((int) disBA*10)/10.0);
			}
			System.out.println(str);
		}

		
		
		
		
		
	} // Of main()
	
	
}
