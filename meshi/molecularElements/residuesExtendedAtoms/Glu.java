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
 *           CG - CD - OE2
 *                | (-)
 *                OE1
 **/
public class Glu extends ResidueExtendedAtoms {
    public final Atom CG, CD, OE1, OE2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       H   CB\n"+
	"           |\n"+
	"           CG - CD - OE2\n"+
	"                |\n"+
	"                OE1\n";
    public Glu(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public Glu(Glu old,int number){
         this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	 setSecondaryStructure(old.secondaryStructure());
	       setAccessibility(old.accessibility());
    }
    public Glu(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(GLU, atomList, number, mode,  addAtomsFlag);

	Object[] temp = new Object[4];
	int i=0;
	temp[i++] = CG = getAtom("CG",ECG, atomList, this);
	temp[i++] = CD = getAtom("CD",ECD, atomList, this);
	temp[i++] = OE1 = getAtom("OE1",EOE, atomList, this);
	temp[i++] = OE2 = getAtom("OE2",EOE, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);

	if ((CB != null)   && (CG != null)) bonds.add(CB.bond(CG));
	if ((CG != null)   && (CD != null)) bonds.add(CG.bond(CD));
	if ((CD != null)   && (OE1 != null)) bonds.add(CD.bond(OE1));
	if ((CD != null)   && (OE2 != null)) bonds.add(CD.bond(OE2));
    }
    public String comment() {
	return COMMENT;
    }
}
