package meshi.parameters;
import meshi.molecularElements.Atom;
public interface AtomTypes {
    //Termini
    public static final int TRN = Atom.addType("TRN"); // Charged Nitrogen at tne N-terminus 
    public static final int TRC = Atom.addType("TRC"); // Carboxyl Carbon at the C-terminus
    public static final int TRO = Atom.addType("TRO"); // Charged Oxygen at tne C-terminus 
    //Ala
    public static final int AH  = Atom.addType("AH");
    public static final int AN  = Atom.addType("AN");
    public static final int ACA = Atom.addType("ACA");
    public static final int AC  = Atom.addType("AC");
    public static final int AO  = Atom.addType("AO");
    public static final int ACB = Atom.addType("ACB");

    //Cys
    public static final int CH  = Atom.addType("CH");
    public static final int CN  = Atom.addType("CN");
    public static final int CCA = Atom.addType("CCA");
    public static final int CC  = Atom.addType("CC");
    public static final int CO  = Atom.addType("CO");
    public static final int CCB = Atom.addType("CCB");
    public static final int CSG = Atom.addType("CSG");

    //Asp
    public static final int DH  = Atom.addType("DH");
    public static final int DN  = Atom.addType("DN");
    public static final int DCA = Atom.addType("DCA");
    public static final int DC  = Atom.addType("DC");
    public static final int DO  = Atom.addType("DO");
    public static final int DCB = Atom.addType("DCB");
    public static final int DCG = Atom.addType("DCG");
    public static final int DOD = Atom.addType("DOD");

    //Glu
    public static final int EH  = Atom.addType("EH");
    public static final int EN  = Atom.addType("EN");
    public static final int ECA = Atom.addType("ECA");
    public static final int EC  = Atom.addType("EC");
    public static final int EO  = Atom.addType("EO");
    public static final int ECB = Atom.addType("ECB");
    public static final int ECG = Atom.addType("ECG");
    public static final int ECD = Atom.addType("ECD");
    public static final int EOE = Atom.addType("EOE");

    //Phe
    public static final int FH  = Atom.addType("FH");
    public static final int FN  = Atom.addType("FN");
    public static final int FCA = Atom.addType("FCA");
    public static final int FC  = Atom.addType("FC");
    public static final int FO  = Atom.addType("FO");
    public static final int FCB = Atom.addType("FCB");
    public static final int FCG = Atom.addType("FCG");
    public static final int FCD = Atom.addType("FCD");
    public static final int FCE = Atom.addType("FCE");
    public static final int FCZ = Atom.addType("FCZ");

    //Gly
    public static final int GH  = Atom.addType("GH");
    public static final int GN  = Atom.addType("GN");
    public static final int GCA = Atom.addType("GCA");
    public static final int GC  = Atom.addType("GC");
    public static final int GO  = Atom.addType("GO");

    //His
    public static final int HH  = Atom.addType("HH");
    public static final int HN  = Atom.addType("HN");
    public static final int HCA = Atom.addType("HCA");
    public static final int HC  = Atom.addType("HC");
    public static final int HO  = Atom.addType("HO");
    public static final int HCB = Atom.addType("HCB");
    public static final int HCG = Atom.addType("HCG");
    public static final int HCD = Atom.addType("HCD");
    public static final int HND = Atom.addType("HND");
    public static final int HHD = Atom.addType("HHD");
    public static final int HCE = Atom.addType("HCE"); 
    public static final int HNE = Atom.addType("HNE");
    public static final int HHE = Atom.addType("HHE");

   //Ile
    public static final int IH  = Atom.addType("IH");
    public static final int IN  = Atom.addType("IN");
    public static final int ICA = Atom.addType("ICA");
    public static final int IC  = Atom.addType("IC");
    public static final int IO  = Atom.addType("IO");
    public static final int ICB = Atom.addType("ICB");
    public static final int ICG1 = Atom.addType("ICG1");
    public static final int ICG2 = Atom.addType("ICG2");
    public static final int ICD = Atom.addType("ICD");

    //Lys
    public static final int KH  = Atom.addType("KH");
    public static final int KN  = Atom.addType("KN");
    public static final int KCA = Atom.addType("KCA");
    public static final int KC  = Atom.addType("KC");
    public static final int KO  = Atom.addType("KO");
    public static final int KCB = Atom.addType("KCB");
    public static final int KCG = Atom.addType("KCG");
    public static final int KCD = Atom.addType("KCD");
    public static final int KCE = Atom.addType("KCE");
    public static final int KNZ = Atom.addType("KNZ");

    //Leu
    public static final int LH  = Atom.addType("LH");
    public static final int LN  = Atom.addType("LN");
    public static final int LCA = Atom.addType("LCA");
    public static final int LC  = Atom.addType("LC");
    public static final int LO  = Atom.addType("LO");
    public static final int LCB = Atom.addType("LCB");
    public static final int LCG = Atom.addType("LCG");
    public static final int LCD1 = Atom.addType("LCD1");
    public static final int LCD2 = Atom.addType("LCD2");

    //Met
    public static final int MH  = Atom.addType("MH");
    public static final int MN  = Atom.addType("MN");
    public static final int MCA = Atom.addType("MCA");
    public static final int MC  = Atom.addType("MC");
    public static final int MO  = Atom.addType("MO");
    public static final int MCB = Atom.addType("MCB");
    public static final int MCG = Atom.addType("MCG");
    public static final int MSD = Atom.addType("MSD");
    public static final int MCE = Atom.addType("MCE");

    //Asn
    public static final int NH  = Atom.addType("NH");
    public static final int NN  = Atom.addType("NN");
    public static final int NCA = Atom.addType("NCA");
    public static final int NC  = Atom.addType("NC");
    public static final int NO  = Atom.addType("NO");
    public static final int NCB = Atom.addType("NCB");
    public static final int NCG = Atom.addType("NCG");
    public static final int NOD = Atom.addType("NOD");
    public static final int NND = Atom.addType("NND");
    public static final int NHD1 = Atom.addType("NHD1");
    public static final int NHD2 = Atom.addType("NHD2");
    //Pro
    public static final int PN  = Atom.addType("PN");
    public static final int PCA = Atom.addType("PCA");
    public static final int PC  = Atom.addType("PC");
    public static final int PO  = Atom.addType("PO");
    public static final int PCB = Atom.addType("PCB");
    public static final int PCG = Atom.addType("PCG");
    public static final int PCD = Atom.addType("PCD");

    //Gln
    public static final int QH  = Atom.addType("QH");
    public static final int QN  = Atom.addType("QN");
    public static final int QCA = Atom.addType("QCA");
    public static final int QC  = Atom.addType("QC");
    public static final int QO  = Atom.addType("QO");
    public static final int QCB = Atom.addType("QCB");
    public static final int QCG = Atom.addType("QCG");
    public static final int QCD = Atom.addType("QCD");
    public static final int QOE = Atom.addType("QOE");
    public static final int QNE = Atom.addType("QNE");
    public static final int QHE1 = Atom.addType("QHE1");
    public static final int QHE2 = Atom.addType("QHE2");

    //Arg
    public static final int RH  = Atom.addType("RH");
    public static final int RN  = Atom.addType("RN");
    public static final int RCA = Atom.addType("RCA");
    public static final int RC  = Atom.addType("RC");
    public static final int RO  = Atom.addType("RO");
    public static final int RCB = Atom.addType("RCB");
    public static final int RCG = Atom.addType("RCG");
    public static final int RCD = Atom.addType("RCD");
    public static final int RNE = Atom.addType("RNE");
    public static final int RHE = Atom.addType("RHE");
    public static final int RCZ = Atom.addType("RCZ");
    public static final int RNH = Atom.addType("RNH");

    //Ser
    public static final int SH  = Atom.addType("SH");
    public static final int SN  = Atom.addType("SN");
    public static final int SCA = Atom.addType("SCA");
    public static final int SC  = Atom.addType("SC");
    public static final int SO  = Atom.addType("SO");
    public static final int SCB = Atom.addType("SCB");
    public static final int SOG = Atom.addType("SOG");

    //Thr
    public static final int TH  = Atom.addType("TH");
    public static final int TN  = Atom.addType("TN");
    public static final int TCA = Atom.addType("TCA");
    public static final int TC  = Atom.addType("TC");
    public static final int TO  = Atom.addType("TO");
    public static final int TCB = Atom.addType("TCB");
    public static final int TCG = Atom.addType("TCG");
    public static final int TOG = Atom.addType("TOG");

    //Val
    public static final int VH  = Atom.addType("VH");
    public static final int VN  = Atom.addType("VN");
    public static final int VCA = Atom.addType("VCA");
    public static final int VC  = Atom.addType("VC");
    public static final int VO  = Atom.addType("VO");
    public static final int VCB = Atom.addType("VCB");
    public static final int VCG1 = Atom.addType("VCG1");
    public static final int VCG2 = Atom.addType("VCG2");

    //Trp
    public static final int WH  = Atom.addType("WH");
    public static final int WN  = Atom.addType("WN");
    public static final int WCA = Atom.addType("WCA");
    public static final int WC  = Atom.addType("WC");
    public static final int WO  = Atom.addType("WO");
    public static final int WCB = Atom.addType("WCB");
    public static final int WCG = Atom.addType("WCG");
    public static final int WCD1 = Atom.addType("WCD1");
    public static final int WCD2 = Atom.addType("WCD2");
    public static final int WCE2 = Atom.addType("WCE2");
    public static final int WCE3 = Atom.addType("WCE3");
    public static final int WNE = Atom.addType("WNE");
    public static final int WHE = Atom.addType("WHE");
    public static final int WCZ2 = Atom.addType("WCZ2");
    public static final int WCZ3 = Atom.addType("WCZ3");
    public static final int WCH2 = Atom.addType("WCH2");
    //Tyr
    public static final int YH  = Atom.addType("YH");
    public static final int YN  = Atom.addType("YN");
    public static final int YCA = Atom.addType("YCA");
    public static final int YC  = Atom.addType("YC");
    public static final int YO  = Atom.addType("YO");
    public static final int YCB = Atom.addType("YCB");
    public static final int YCG = Atom.addType("YCG");
    public static final int YCD = Atom.addType("YCD");
    public static final int YCE = Atom.addType("YCE");
    public static final int YCZ = Atom.addType("YCZ");
    public static final int YOH = Atom.addType("YOH","DONE");
    //
    public static final int[] BB_HYDROGENS = {AH, CH, DH, EH, FH, GH, HH, IH, KH, LH, MH, NH, QH, RH, SH, TH, VH, WH, YH};
    public static final int[] BB_OXYGENS = {AO, CO, DO, EO, FO, GO, HO, IO, KO, LO, MO, NO, PO, QO, RO, SO, TO, VO, WO, YO};    
    public static final int[] BB_NITROGENS = {AN, CN, DN, EN, FN, GN, HN, IN, KN, LN, MN, NN, PN, QN, RN, SN, TN, VN, WN,YN, PN, TRN};
    public static final int[] BB_CARBONS = {AC, CC, DC, EC, FC, GC, HC, IC, KC, LC, MC, NC, PC, QC, RC, SC, TC, VC, WC,YC, PC, TRC};} 
    
