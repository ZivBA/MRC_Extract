package meshi.util.crossLinking;

import meshi.molecularElements.AtomList;
import meshi.util.TRIC.PutUnitInAnyTopPositionYeast;

public class CrosslinkVectorPIC extends CrosslinkVector {

	private static final long serialVersionUID = 1L;

	public CrosslinkVectorPIC() {
		// TODO Auto-generated constructor stub
	}

	public CrosslinkVectorPIC(String filename , int fileType) {
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
		return new CrosslinkVectorPIC();
	}	
	
	public static void main(String[] args) {


		CrosslinkVectorPIC xlScience = new CrosslinkVectorPIC("H:\\PIC_new_map\\ModelForDeposit\\Science_H_dataset.txt",4);
		AtomList PICmodel = new AtomList("H:\\PIC_new_map\\ModelForDeposit\\Form1\\Murakami_PIC_Form1_Alt_Ssl1_pos.pdb");
		System.out.println(xlScience.printDistancesOnStructure(PICmodel, 10));
		System.exit(0);		
		
//		CrosslinkVectorPIC xlVec11 = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\PIC_SUB1_Jeff\\MATLAB_Trypsin\\PIC_JL112013_120913_Trypsin_pls_SUB1_All_Pooled_High_confidence_Paper_Format.txt",4);
//		CrosslinkVectorPIC xlVec22 = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\PIC_SUB1_Jeff\\MATLAB_Trypsin\\Jeff_in_paper_format_with_distances.txt",4);
//		System.out.print(xlVec11.compareToOtherXLvec(xlVec22));
//		System.exit(0);
				
//		CrosslinkVectorPIC xlVec11 = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\PIC_03_2014\\Science_S_data_for_comparisons.txt",4);
//		CrosslinkVectorPIC xlVec22 = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\PIC_03_2014\\Science_K_data_for_comparisons.txt",4);
//		xlVec11.addVecNonRedundant(xlVec22);
//		CrosslinkVectorPIC xlVecNew = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\PIC_03_2014\\MATLAB_files_trypsin_HCD\\PIC_140319_trypsin_HCD_SUB1_all_fraction_pooled_just_High_Confidence.txt",4);
//		System.out.print(xlVec11.compareToOtherXLvec(xlVecNew));
//		System.exit(0);
				
//		CrosslinkVectorPIC xlVec = new CrosslinkVectorPIC("G:\\PIC_SUB1\\Jeff_in_paper_format.txt",4);
//		CrosslinkVectorPIC xlVec1 = (CrosslinkVectorPIC) (new CrosslinkVectorPIC("G:\\PIC_SUB1\\Science_S_K_sets.txt",4)).filterOutRedundancies();
//		xlVec1.setCreatingFile("Science_K_and_S");
//		System.out.print(xlVec1.compareToOtherXLvec(xlVec));
//		//System.out.print(xlVec.printDistancesOnStructure(new AtomList("C:\\Users\\Nir\\PIC\\HomologyModeling\\AllTogether\\All_together.pdb"), 5));
//		System.exit(0);

		CrosslinkVectorPIC xlVecKq = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\PIC_03_2013\\PIC_S_trypsin\\CLOSED_S_2012_with_distances.txt",1);
		CrosslinkVectorPIC xlVecMergedSq = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\PIC_03_2013\\PIC_S_trypsin\\CLOSED_S_Trypsin_03_2013_with_distances_cut_FPR_15.txt",1);
		System.out.print(xlVecKq.compareToOtherXLvec(xlVecMergedSq));
		System.exit(0);
		

		
		CrosslinkVectorPIC xlVecSdup = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\PIC_03_2013\\Files_for_Science_revision\\CLOSED_S_PIC_XL_table_merged_duplicate.txt",1);
		System.out.println(xlVecSdup.filterOutRedundancies().size());
		System.out.println(xlVecSdup.filterOutRedundancies().filterOutInterUnit().size());
		System.out.println(xlVecSdup.filterOutRedundancies().filterOutIntraUnit().size());
		CrosslinkVectorPIC xlVecKdup = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\PIC_03_2013\\Files_for_Science_revision\\PAPER_K_PIC_table.txt",1);	
		xlVecSdup.addVec(xlVecKdup);
		System.out.println(xlVecSdup.filterOutRedundancies().size());
		for (Crosslink xl : xlVecSdup.filterOutRedundancies()) {
			if ( ((xl.protName1().equals("Tfg1") | xl.protName1().equals("Tfg2")) & !(xl.protName2().equals("Tfg1") | xl.protName2().equals("Tfg2"))) |
					(!(xl.protName1().equals("Tfg1") | xl.protName1().equals("Tfg2")) & (xl.protName2().equals("Tfg1") | xl.protName2().equals("Tfg2")))   ){
				System.out.print(xl);
			}
		}

				
		System.exit(0);
		
		CrosslinkVectorPIC xlVecK = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\PIC_03_2013\\PIC_S_trypsin\\CLOSED_K_for_comparison_with_merged_S.txt",1);
		CrosslinkVectorPIC xlVecMergedS = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\PIC_03_2013\\PIC_S_trypsin\\CLOSED_S_2012_and_2013_with_distances.txt",1);
		for (Crosslink xl : xlVecMergedS) {
			if (xlVecK.find(xl.protName1(), xl.absPos1(), xl.protName2(), xl.absPos2()) != null) {
				System.out.println("Yes");
			}
			else {
				if (xl.protName1().equals("1") | xl.protName1().equals("3") | xl.protName1().equals("4") |
						xl.protName2().equals("1") | xl.protName2().equals("3") | xl.protName2().equals("4")) {
					System.out.println("N/A");
				}
				else {
					System.out.println("No");
				}
			}
		}
		
//		CrosslinkVectorPIC xlVec = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\PIC_03_2013\\PIC_S_trypsin\\CLOSED_S_Trypsin_03_2013.txt",1);
//		System.out.print(xlVec.printDistancesOnStructure(new AtomList("C:\\Users\\Nir\\PIC\\HomologyModeling\\AllTogether\\All_together.pdb"), 5));
		
		
		CrosslinkVectorPIC xlVecSold = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\PIC_03_2013\\PIC_S_trypsin\\Paper_CLOSED_S_FPR15.txt",1);
		CrosslinkVectorPIC xlVecSnew = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\PIC_03_2013\\PIC_S_trypsin\\CLOSED_S_03_2013_FPR15.txt",1);
		//System.out.print(xlVecSold.compareToOtherXLvec(xlVecSnew));
		
		System.out.println("In old:\n----------------");
		for (Crosslink xl : xlVecSold) {
			if (xl.protName1().equals("1") && xl.protName2().equals("2")) {
				System.out.println(xl);
			}
		}
//		System.out.println(xlVecSold.filterProteinsInSet("F"));
		System.out.println("\n\n\nIn new:\n----------------");
		for (Crosslink xl : xlVecSnew) {
			if (xl.protName1().equals("1") && xl.protName2().equals("2")) {
				System.out.println(xl);
			}
		}
//		System.out.println(xlVecSnew.filterProteinsInSet("F"));
		

		CrosslinkVectorPIC xlVecS = new CrosslinkVectorPIC("C:\\Users\\Nir\\myIMP\\automatic_combinatorial\\S_dataset.txt",1);
		System.out.println(xlVecS.size());
//		CrosslinkVectorPIC xlVecK = new CrosslinkVectorPIC("C:\\Users\\Nir\\myIMP\\automatic_combinatorial\\K_dataset.txt",1);
		System.out.println(xlVecK.size());
		xlVecS.addVec(xlVecK);
		System.out.println(xlVecS.filterOutRedundancies().size());
		System.out.println(xlVecS.filterOutRedundancies().filterOutInterUnit().size());
		System.out.println(xlVecS.filterOutRedundancies().filterOutIntraUnit().size());
		System.exit(0);
		System.out.println(xlVecS.filterOutRedundancies().filterOutProteinsInSet("ABCDEFGHIJKLMNOPSTU6").filterOutIntraUnit().size());

//		CrosslinkVectorPIC xlVecSold = new CrosslinkVectorPIC("C:\\Users\\Nir\\myIMP\\automatic_combinatorial\\S_dataset.txt",1);
//		CrosslinkVectorPIC xlVecSnew = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\PIC_01_2013\\PIC_S_Trypsin_inclu\\PIC_S_Trypsin_Inclu_List_all_fractions_pooled_high_confidence.txt",1);
//		System.out.print(xlVecSold.compareToOtherXLvec(xlVecSnew));
//		System.exit(0);
		
				
//		System.out.println(xlVecS.compareToOtherXLvec(xlVecK));
		System.out.print(xlVecS.filterProteinsInSet("2").filterOutInterUnit());
		System.out.print(xlVecK.filterProteinsInSet("2").filterOutInterUnit());

		
		
//		CrosslinkVectorPIC xlVec = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\AnalysisCenter\\PIC_H_FOR_PAPER\\H_120511_ALL_POOLED.txt",1);
//		System.out.print(xlVec.printDistancesOnStructure(new AtomList("C:\\Users\\Nir\\PIC\\HomologyModeling\\AllTogether\\All_together.pdb"), 5));

//		CrosslinkVectorPIC xlVec = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\AnalysisCenter\\PIC_CLOSED_K_FOR_PAPER\\CLOSED_K_ALL_POOLED.txt",1);
//		System.out.print(xlVec.printDistancesOnStructure(new AtomList("C:\\Users\\Nir\\PIC\\HomologyModeling\\AllTogether\\All_together.pdb"), 5));

//		CrosslinkVectorPIC xlVec = new CrosslinkVectorPIC("C:\\Users\\Nir\\PIC\\AnalysisCenter\\PIC_CLOSED_S_FOR_PAPER\\CLOSED_S_ALL_POOLED.txt",1);
//		System.out.print(xlVec.printDistancesOnStructure(new AtomList("C:\\Users\\Nir\\PIC\\HomologyModeling\\AllTogether\\All_together.pdb"), 5));
		
//		CrosslinkVectorPIC xlVec = new CrosslinkVectorPIC("tmp_S_set_medium_confidence.txt",1);
//		String gene2letter = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456";
//		int[][] conn = new int[32][32];
//		for (int c=0 ; c<gene2letter.length() ; c++) {
//			for (int d=0 ; d<gene2letter.length() ; d++) {
//				for (Crosslink xl : xlVec) {
//					if ((xl.protName1().equals(gene2letter.substring(c, c+1)) && xl.protName2().equals(gene2letter.substring(d, d+1))) ||
//							(xl.protName2().equals(gene2letter.substring(c, c+1)) && xl.protName1().equals(gene2letter.substring(d, d+1)))) {
//						conn[c][d]++;
//					}
//				}
//			}
//		}
//		System.out.print(" | ");
//		for (int c=0 ; c<gene2letter.length() ; c++) {
//			System.out.print(gene2letter.substring(c, c+1) + "  ");
//		}
//		System.out.println("\n-------------------------------------------------------------------------------------------------------------------------------------------------");
//		for (int c=0 ; c<gene2letter.length() ; c++) {
//			System.out.print(gene2letter.substring(c, c+1) + "|");
//			for (int d=0 ; d<gene2letter.length() ; d++) {
//				if (conn[c][d]<10) {
//					System.out.print(" " + conn[c][d] + " ");
//				}
//				else {
//					System.out.print(conn[c][d] + " ");
//				}
//			}
//			System.out.println();
//		}
		
		
		
		//System.out.println(xlVec1.compareToOtherXLvec(xlVec2));		
	} // Of main()
	
	
}
