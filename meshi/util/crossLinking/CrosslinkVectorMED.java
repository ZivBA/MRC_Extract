package meshi.util.crossLinking;

import meshi.molecularElements.AtomList;

public class CrosslinkVectorMED extends CrosslinkVector {

	public CrosslinkVectorMED() {
		// TODO Auto-generated constructor stub
	}

	public CrosslinkVectorMED(String filename , int fileType) {
		super(filename ,  fileType);
		// TODO Auto-generated constructor stub
	}

	/**
	 * In PIC the stochiometry of all the subunits is 1.
	 */
	public String identicalChains(String protName) {
		return protName;
	}

	public CrosslinkVector createEmptyVector() {
		return new CrosslinkVectorMED();
	}	
	
	public static void main(String[] args) {

		CrosslinkVectorMED xlVec = new CrosslinkVectorMED("C:\\Users\\Nir\\Mediator\\holo_xlink_files\\holo_xlink_files\\xlink_datasets\\bs3_sec_dataset_forMESHI.txt",3);
		System.out.println(xlVec.printDistancesOnStructure(new AtomList("C:\\Users\\Nir\\Mediator\\holo_xlink_files\\holo_xlink_files\\PDB_structures\\1WCM_12subunits.pdb"), 20));
		
		System.out.println(xlVec.size() + " " + xlVec.filterOutInterUnit().size() + " " + xlVec.filterOutIntraUnit().size());
		//System.out.println(xlVec.filterOutBothProteinsInSet("ABCDEFGHIJKL").filterOutBothProteinsInSet("MNOPQRSTUVWXYZ123456"));
		System.out.println(xlVec.filterProteinsInSet("PX").filterOutIntraUnit());
		
	} // Of main()
	
	
}
