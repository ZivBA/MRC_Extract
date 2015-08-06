package meshi.molecularElements.residuesExtendedAtoms;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
/**
 *<pre>
 *                O
 *                |
 *       N - CA - C...n
 *       |   |
 *       H   CB             
 *           |              
 *           CG - CD - CE - NZ
 *
 **/
public class Lys extends ResidueExtendedAtoms {
    public final Atom CG, CD, CE, NZ;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       H   CB\n"+             
	"           |\n"+              
	"           CG - CD - CE - NZ\n";
    public Lys(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public Lys(Lys old,int number){
         this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	 setSecondaryStructure(old.secondaryStructure());
	       setAccessibility(old.accessibility());
    }
    public Lys(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(LYS, atomList, number, mode,  addAtomsFlag);

	Object[] temp = new Object[4];
	int i=0;
	temp[i++] = CG = getAtom("CG",KCG, atomList, this);
	temp[i++] = CD = getAtom("CD",KCD, atomList, this);
	temp[i++] = CE = getAtom("CE",KCE, atomList, this);
	temp[i++] = NZ = getAtom("NZ",KNZ, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);

	if ((CB  != null)   && (CG  != null)) bonds.add(CB.bond(CG));
	if ((CG  != null)   && (CD  != null)) bonds.add(CG.bond(CD));
	if ((CD  != null)   && (CE  != null)) bonds.add(CD.bond(CE));
	if ((CE  != null)   && (NZ  != null)) bonds.add(CE.bond(NZ));
    }
    public String comment() {
	return COMMENT;
    }
}
