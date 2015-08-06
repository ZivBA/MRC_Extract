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
 *           CG - CD - NE2 - HE21
 *                |    |
 *                OE1  HE22
 **/
public class Gln extends ResidueExtendedAtoms {
    public final Atom CG, CD, OE1, NE2, HE21, HE22;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       H   CB\n"+
	"           |\n"+
	"           CG - CD - NE2 - HE21\n"+
	"                |    |\n"+
	"                OE1  HE22\n";
    public Gln(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public Gln(Gln old,int number){
         this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	 setSecondaryStructure(old.secondaryStructure());
	       setAccessibility(old.accessibility());
       }
    public Gln(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(GLN, atomList, number, mode,  addAtomsFlag);
	Object[] temp = new Object[6];
	int i=0;
	temp[i++] = CG = getAtom("CG",QCG, atomList, this);
	temp[i++] = CD = getAtom("CD",QCD, atomList, this);
	temp[i++] = NE2 = getAtom("NE2",QNE, atomList, this);
	if (NE2 != null) {
	    temp[i++] = HE21 = getAtom("HE21",QHE1, atomList, this);
	    temp[i++] = HE22 = getAtom("HE22",QHE2, atomList, this);
	}
	else HE21 = HE22 = null;
	    
	temp[i++] = OE1 = getAtom("OE1",QOE, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);

	if ((CB != null)   && (CG != null))bonds.add(CB.bond(CG));
	if ((CG != null)   && (CD != null))bonds.add(CG.bond(CD));
	if ((CD != null)   && (OE1 != null))bonds.add(CD.bond(OE1));
	if ((CD != null)   && (NE2 != null))bonds.add(CD.bond(NE2));
	if ((NE2 != null)   && (HE21 != null))bonds.add(NE2.bond(HE21));
	if ((NE2 != null)   && (HE22 != null))bonds.add(NE2.bond(HE22));
     }
    public String comment() {
	return COMMENT;
    }
}
