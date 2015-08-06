package meshi.molecularElements.residuesExtendedAtoms;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
/**
 *<pre>
 *                O
 *                |
 *       N - CA - C...n
 *       |   |
 *       H   CB - CG - CD1
 *                |
 *                CD2
 **/
public class Leu extends ResidueExtendedAtoms {
    public final Atom CG, CD1, CD2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       H   CB - CG - CD1\n"+
	"                |\n"+
	"                CD2\n";
    public Leu(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public Leu(Leu old,int number){
         this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	 setSecondaryStructure(old.secondaryStructure());
	       setAccessibility(old.accessibility());
    }
    public Leu(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(LEU, atomList, number, mode,  addAtomsFlag);

	Object[] temp = new Object[3];
	int i=0;
	temp[i++] = CG = getAtom("CG",LCG, atomList, this);
	temp[i++] = CD1 = getAtom("CD1",LCD1, atomList, this);
	temp[i++] = CD2 = getAtom("CD2",LCD2, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);

	if ((CB  != null)   && (CG  != null)) bonds.add(CB.bond(CG));
	if ((CG  != null)   && (CD1  != null)) bonds.add(CG.bond(CD1));
	if ((CG  != null)   && (CD2  != null)) bonds.add(CG.bond(CD2));
    }
    public String comment() {
	return COMMENT;
    }
}
