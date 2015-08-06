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
 *           CG - CD1 - CE1
 *           |          |
 *           CD2- CE2 - CZ-OH
 **/
public class Tyr extends ResidueExtendedAtoms {
    public final Atom CG, CD1, CE1, CZ, CD2, CE2, OH;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       H   CB\n"+
	"           |\n"+
	"           CG - CD1 - CE1\n"+
	"           |          |\n"+
	"           CD2- CE2 - CZ-OH\n";
    public Tyr(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public Tyr(Tyr old,int number){
	this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	setSecondaryStructure(old.secondaryStructure());
	      setAccessibility(old.accessibility());
    }
   public Tyr(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(TYR, atomList, number, mode,  addAtomsFlag);

	Object[] temp = new Object[7];
	int i=0;
	temp[i++] = CG = getAtom("CG",YCG, atomList, this);
	temp[i++] = CD1 = getAtom("CD1",YCD, atomList, this);
	temp[i++] = CE1 = getAtom("CE1",YCE, atomList, this);
	temp[i++] = CZ = getAtom("CZ",YCZ, atomList, this);
	temp[i++] = CD2 = getAtom("CD2",YCD, atomList, this);
	temp[i++] = CE2 = getAtom("CE2",YCE, atomList, this);
	temp[i++] = OH = getAtom("OH",YOH, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);

	if ((CB  != null)  && (CG  != null)) bonds.add(CB.bond(CG));
	if ((CG  != null)  && (CD1 != null)) bonds.add(CG.bond(CD1));
	if ((CG  != null)  && (CD2 != null)) bonds.add(CG.bond(CD2));
	if ((CD1 != null)  && (CE1 != null)) bonds.add(CD1.bond(CE1));
	if ((CD2 != null)  && (CE2 != null)) bonds.add(CD2.bond(CE2));
	if ((CE1 != null)  && (CZ  != null)) bonds.add(CE1.bond(CZ));
	if ((CE2 != null)  && (CZ  != null)) bonds.add(CE2.bond(CZ));
	if ((CZ  != null)  && (OH  != null)) bonds.add(CZ.bond(OH));
    }
    public String comment() {
	return COMMENT;
    }
}

