package meshi.molecularElements.residuesExtendedAtoms;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
/**
 *<pre>
 *                O
 *                |
 *       N - CA - C...n
 *       |   |
 *       H   CB   CE3 - CZ3 - CH2
 *           |    |           |
 *           CG - CD2 - CE2 - CZ2
 *           |           |
 *           CD1 ------ NE1 - HE1
 *
 **/
public class Trp extends ResidueExtendedAtoms {
    public final Atom CG, CD1, NE1, HE1, CE2, CZ2, CH2, CZ3, CE3, CD2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       H   CB   CE3 - CZ3 - CH2\n"+
	"           |    |           |\n"+
	"           CG - CD2 - CE2 - CZ2\n"+
	"           |           |\n"+
	"           CD1 ------ NE1 - HE1\n";
    public Trp(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public Trp(Trp old,int number){
         this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	 setSecondaryStructure(old.secondaryStructure());
	       setAccessibility(old.accessibility());
    }
   public Trp(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(TRP, atomList, number, mode,  addAtomsFlag);
	Object[] temp = new Object[10];
	int i=0;
	temp[i++] = CG = getAtom("CG",WCG, atomList, this);
	temp[i++] = CD1 = getAtom("CD1",WCD1, atomList, this);
	temp[i++] = NE1 = getAtom("NE1",WNE, atomList, this);
	if (NE1 != null)
	    temp[i++] = HE1 = getAtom("HE1",WHE, atomList, this);
	else HE1 = null;
	temp[i++] = CD2 = getAtom("CD2",WCD2, atomList, this);

	temp[i++] = CE2 = getAtom("CE2",WCE2, atomList, this);
	temp[i++] = CE3 = getAtom("CE3",WCE3, atomList, this);
	temp[i++] = CZ3 = getAtom("CZ3",WCZ3, atomList, this);
	temp[i++] = CH2 = getAtom("CH2",WCH2, atomList, this);
	temp[i++] = CZ2 = getAtom("CZ2",WCZ2, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);

	if ((CB  != null)   && (CG   != null)) bonds.add(CB.bond(CG));
 	if ((CG  != null)   && (CD1  != null)) bonds.add(CG.bond(CD1));
 	if ((CG  != null)   && (CD2  != null)) bonds.add(CG.bond(CD2));
 	if ((CD1 != null)   && (NE1  != null)) bonds.add(CD1.bond(NE1));
 	if ((CD2 != null)   && (CE2  != null)) bonds.add(CD2.bond(CE2));
 	if ((NE1 != null)   && (HE1  != null)) bonds.add(NE1.bond(HE1));
 	if ((CE2 != null)   && (CZ2  != null)) bonds.add(CE2.bond(CZ2));
 	if ((CZ2 != null)   && (CH2  != null)) bonds.add(CZ2.bond(CH2));
 	if ((CH2 != null)   && (CZ3  != null)) bonds.add(CH2.bond(CZ3));
 	if ((CZ3 != null)   && (CE3  != null)) bonds.add(CZ3.bond(CE3));
 	if ((CE3 != null)   && (CD2  != null)) bonds.add(CE3.bond(CD2));
 	if ((NE1 != null)   && (CE2  != null)) bonds.add(NE1.bond(CE2));
      }
    public String comment() {
	return COMMENT;
    }
}
