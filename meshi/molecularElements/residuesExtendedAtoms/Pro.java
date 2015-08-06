package meshi.molecularElements.residuesExtendedAtoms;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
/**
 *<pre>
 *                O
 *                |
 *       N - CA - C...n
 *       |   |
 *       CD  CB
 *       \  /
 *        CG
 **/
public class Pro extends ResidueExtendedAtoms {
    public final Atom CG, CD;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       CD  CB\n"+
	"       \\  /\n"+
	"        CG\n";
    public Pro(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public Pro(Pro old,int number){
         this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	 setSecondaryStructure(old.secondaryStructure());
	       setAccessibility(old.accessibility());
    }
    public Pro(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(PRO, atomList, number, mode,  addAtomsFlag);

	Object[] temp = new Object[2];
	int i=0;
	temp[i++] = CG = getAtom("CG",PCG, atomList, this);
	temp[i++] = CD = getAtom("CD",PCD, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);

	if ((CB  != null)   && (CG  != null)) bonds.add(CB.bond(CG));
	if ((CG  != null)   && (CD  != null)) bonds.add(CG.bond(CD));
	if ((CD  != null)   && (N  != null)) bonds.add(CD.bond(N));
    }
    public String comment() {
	return COMMENT;
    }
}
