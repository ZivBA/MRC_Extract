package meshi.parameters;
import meshi.molecularElements.Residue;
public interface Residues extends AtomTypes{
    //Termini
    public static final int ALA = Residue.addName("ALA","A",ACA); 
    public static final int CYS = Residue.addName("CYS","C",CCA); 
    public static final int ASP = Residue.addName("ASP","D",DCA); 
    public static final int GLU = Residue.addName("GLU","E",ECA); 
    public static final int PHE = Residue.addName("PHE","F",FCA); 
    public static final int GLY = Residue.addName("GLY","G",GCA); 
    public static final int HIS = Residue.addName("HIS","H",HCA); 
    public static final int ILE = Residue.addName("ILE","I",ICA); 
    public static final int LYS = Residue.addName("LYS","K",KCA); 
    public static final int LEU = Residue.addName("LEU","L",LCA); 
    public static final int MET = Residue.addName("MET","M",MCA); 
    public static final int ASN = Residue.addName("ASN","N",NCA); 
    public static final int PRO = Residue.addName("PRO","P",PCA); 
    public static final int GLN = Residue.addName("GLN","Q",QCA); 
    public static final int ARG = Residue.addName("ARG","R",RCA); 
    public static final int SER = Residue.addName("SER","S",SCA); 
    public static final int THR = Residue.addName("THR","T",TCA); 
    public static final int VAL = Residue.addName("VAL","V",VCA); 
    public static final int TRP = Residue.addName("TRP","W",WCA); 
    public static final int TYR = Residue.addName("TYR","Y",YCA);
    public static final int UNK = Residue.addName("UNK","X",-1,"DONE");
    public static final int SINGLE = 0;
    public static final int NTER = 1;
    public static final int NORMAL = 2;
    public static final int CTER = 3;
    public static final int ADD_ATOMS_AND_FREEZE = 0;
    public static final int ADD_BACKBONE_AND_FREEZE = 4;
    public static final int ADD_ATOMS = 3;
    public static final int DO_NOT_ADD_ATOMS = 1;
    public static final int ADD_HYDROGENS_AND_FREEZE = 2;
    public static final int ADD_NOT_SET = -1;
}
 
  
