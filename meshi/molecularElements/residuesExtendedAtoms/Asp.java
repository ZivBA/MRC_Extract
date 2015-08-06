package meshi.molecularElements.residuesExtendedAtoms;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
/**
 *<pre>
 *
 *                O
 *                |
 *       N - CA - C...n
 *       |   |
 *       H   CB
 *           |
 *           CG - OD2
 *           | (-)
 *           OD1
 **/
public class Asp extends ResidueExtendedAtoms {
    public final Atom CG, OD1, OD2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       H   CB\n"+
	"           |\n"+
	"           CG - OD2\n"+
	"           | \n"+
	"           OD1\n";
    public Asp(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public Asp(Asp old,int number){
         this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	 setSecondaryStructure(old.secondaryStructure());
	       setAccessibility(old.accessibility());
       }
    public Asp(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(ASP, atomList, number, mode,  addAtomsFlag);

	Object[] temp = new Object[3];
	int i=0;
	temp[i++] = CG = getAtom("CG",DCG, atomList, this);
	temp[i++] = OD1 = getAtom("OD1",DOD, atomList, this);
	temp[i++] = OD2 = getAtom("OD2",DOD, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);

	if ((CB != null)   && (CG != null))  bonds.add(CB.bond(CG));
	if ((CG != null)   && (OD2 != null))  bonds.add(CG.bond(OD2));
	if ((CG != null)   && (OD1 != null))  bonds.add(CG.bond(OD1));
    }
    public String comment() {
	return COMMENT;
    }
}
