package meshi.molecularElements;

public class MESHI2TsaiConversion {
	/** Converting from MESHI's 190 atom types to Tsai's 14. This is according to Tsai 99'.
	 * The Tsai types are:
	 * 0 - Nitrogen with zero hydrogens bonded (only PN)
	 * 1 - Nitrogen with 1 hydrogens bonded (e.g. AN)
	 * 2 - Nitrogen with 2 hydrogens bonded (e.g. NND)
	 * 3 - Nitrogen with 3 hydrogens bonded (e.g. KNZ)
	 * 4 - Oxygen with 0 hydrogens bonded - charged (e.g. DOD)
	 * 5 - Oxygen with 1 hydrogens bonded (e.g. SOG)
	 * 6 - Carbon with 0 hydrogens bonded (e.g. AC)
	 * 7 - Carbon with 1 hydrogens bonded aromatic (e.g. FCD)
	 * 8 - Carbon with 1 hydrogens bonded but not aromatic (e.g. ACA)
	 * 9 - Carbon with 2 hydrogens bonded (e.g. CCB)
	 * 10 - Carbon with 3 hydrogens bonded (e.g. ACB)
	 * 11 - MSD Sulfor 
	 * 12 - CSG Sulfor 
	 * 13 - any Hydrogen
	 */
	public static final int[] MESHI2Tsai = {
		3,   /* TRN  */
		6,   /* TRC  */
		4,   /* TRO  */
		13,   /* AH   */
		1,   /* AN   */
		8,   /* ACA  */
		6,   /* AC   */
		4,   /* AO   */
		10,   /* ACB  */
		13,   /* CH   */
		1,   /* CN   */
		8,   /* CCA  */
		6,   /* CC   */
		4,   /* CO   */
		9,   /* CCB  */
		12,   /* CSG  */
		13,   /* DH   */
		1,   /* DN   */
		8,   /* DCA  */
		6,   /* DC   */
		4,   /* DO   */
		9,   /* DCB  */
		6,   /* DCG  */
		5,   /* DOD  */
		13,   /* EH   */
		1,   /* EN   */
		8,   /* ECA  */
		6,   /* EC   */
		4,   /* EO   */
		9,   /* ECB  */
		9,   /* ECG  */
		6,   /* ECD  */
		5,   /* EOE  */
		13,   /* FH   */
		1,   /* FN   */
		8,   /* FCA  */
		6,   /* FC   */
		4,   /* FO   */
		9,   /* FCB  */
		6,   /* FCG  */
		7,   /* FCD  */
		7,   /* FCE  */
		7,   /* FCZ  */
		13,   /* GH   */
		1,   /* GN   */
		9,   /* GCA  */
		6,   /* GC   */
		4,   /* GO   */
		13,   /* HH   */
		1,   /* HN   */
		8,   /* HCA  */
		6,   /* HC   */
		4,   /* HO   */
		9,   /* HCB  */
		6,   /* HCG  */
		7,   /* HCD  */
		1,   /* HND  */
		13,   /* HHD  */
		7,   /* HCE  */
		1,   /* HNE  */
		13,   /* HHE  */
		13,   /* IH   */
		1,   /* IN   */
		8,   /* ICA  */
		6,   /* IC   */
		4,   /* IO   */
		8,   /* ICB  */
		9,   /* ICG1 */
		10,   /* ICG2 */
		10,   /* ICD  */
		13,   /* KH   */
		1,   /* KN   */
		8,   /* KCA  */
		6,   /* KC   */
		4,   /* KO   */
		9,   /* KCB  */
		9,   /* KCG  */
		9,   /* KCD  */
		9,   /* KCE  */
		3,   /* KNZ  */
		13,   /* LH   */
		1,   /* LN   */
		8,   /* LCA  */
		6,   /* LC   */
		4,   /* LO   */
		9,   /* LCB  */
		8,   /* LCG  */
		10,   /* LCD1 */
		10,   /* LCD2 */
		13,   /* MH   */
		1,   /* MN   */
		8,   /* MCA  */
		6,   /* MC   */
		4,   /* MO   */
		9,   /* MCB  */
		9,   /* MCG  */
		11,   /* MSD  */
		10,   /* MCE  */
		13,   /* NH   */
		1,   /* NN   */
		8,   /* NCA  */
		6,   /* NC   */
		4,   /* NO   */
		9,   /* NCB  */
		6,   /* NCG  */
		4,   /* NOD  */
		2,   /* NND  */
		13,   /* NHD1 */
		13,   /* NHD2 */
		0,   /* PN   */
		8,   /* PCA  */
		6,   /* PC   */
		4,   /* PO   */
		9,   /* PCB  */
		9,   /* PCG  */
		9,   /* PCD  */
		13,   /* QH   */
		1,   /* QN   */
		8,   /* QCA  */
		6,   /* QC   */
		4,   /* QO   */
		9,   /* QCB  */
		9,   /* QCG  */
		6,   /* QCD  */
		4,   /* QOE  */
		2,   /* QNE  */
		13,   /* QHE1 */
		13,   /* QHE2 */
		13,   /* RH   */
		1,   /* RN   */
		8,   /* RCA  */
		6,   /* RC   */
		4,   /* RO   */
		9,   /* RCB  */
		9,   /* RCG  */
		9,   /* RCD  */
		1,   /* RNE  */
		13,   /* RHE  */
		6,   /* RCZ  */
		3,   /* RNH  */
		13,   /* SH   */
		1,   /* SN   */
		8,   /* SCA  */
		6,   /* SC   */
		4,   /* SO   */
		9,   /* SCB  */
		5,   /* SOG  */
		13,   /* TH   */
		1,   /* TN   */
		8,   /* TCA  */
		6,   /* TC   */
		4,   /* TO   */
		8,   /* TCB  */
		10,   /* TCG  */
		5,   /* TOG  */
		13,   /* VH   */
		1,   /* VN   */
		8,   /* VCA  */
		6,   /* VC   */
		4,   /* VO   */
		8,   /* VCB  */
		10,   /* VCG1 */
		10,   /* VCG2 */
		13,   /* WH   */
		1,   /* WN   */
		8,   /* WCA  */
		6,   /* WC   */
		4,   /* WO   */
		9,   /* WCB  */
		6,   /* WCG  */
		7,   /* WCD1 */
		6,   /* WCD2 */
		6,   /* WCE2 */
		7,   /* WCE3 */
		1,   /* WNE  */
		13,   /* WHE  */
		7,   /* WCZ2 */
		7,   /* WCZ3 */
		7,   /* WCH2 */
		13,   /* YH   */
		1,   /* YN   */
		8,   /* YCA  */
		6,   /* YC   */
		4,   /* YO   */
		9,   /* YCB  */
		6,   /* YCG  */
		7,   /* YCD  */
		7,   /* YCE  */
		6,   /* YCZ  */
		5   /* YOH  */	};		
	
}
