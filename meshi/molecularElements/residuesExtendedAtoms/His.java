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
 *           CG - CD2 - NE2 - HE2
 *           |           |
 *     HD1 - ND1 ------ CE1
 **/
public class His extends  ResidueExtendedAtoms{
    public final Atom CG, CD2, NE2, HE2, CE1, ND1, HD1;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"               O\n"+
	"               |\n"+
	"      N - CA - C...n\n"+
	"      |   |\n"+
	"      H   CB\n"+
	"          |\n"+
	"          CG - CD2 - NE2 - HE2\n"+
	"          |           |\n"+
	"    HD1 - ND1 ------ CE1\n";
    public His(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public His(His old,int number){
         this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	 setSecondaryStructure(old.secondaryStructure());
	       setAccessibility(old.accessibility());
    }
    public His(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(HIS, atomList, number, mode,  addAtomsFlag);

	Object[] temp = new Object[7];
	int i=0;
	temp[i++] = CG = getAtom("CG",HCG, atomList, this);
	temp[i++] = ND1 = getAtom("ND1",HND, atomList, this);
	if (ND1 != null) 
	    temp[i++] = HD1 = getAtom("HD1",HHD, atomList, this);
	else HD1 = null;
	temp[i++] = CE1 = getAtom("CE1",HCE, atomList, this);
	temp[i++] = NE2 = getAtom("NE2",HNE, atomList, this);
	if (NE2 != null) 
	    temp[i++] = HE2 = getAtom("HE2",HHE, atomList, this);
	else HE2 = null;
	temp[i++] = CD2 = getAtom("CD2",HCD, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);

	if ((CB  != null)   && (CG  != null)) bonds.add(CB.bond(CG));
 	if ((CG  != null)   && (ND1 != null)) bonds.add(CG.bond(ND1));
 	if ((ND1 != null)   && (HD1 != null)) bonds.add(ND1.bond(HD1));
 	if ((ND1 != null)   && (CE1 != null)) bonds.add(ND1.bond(CE1));
 	if ((CG  != null)   && (CD2 != null)) bonds.add(CG.bond(CD2));
 	if ((CD2 != null)   && (NE2 != null)) bonds.add(CD2.bond(NE2));
 	if ((CE1 != null)   && (NE2 != null)) bonds.add(CE1.bond(NE2));
 	if ((NE2 != null)   && (HE2 != null)) bonds.add(NE2.bond(HE2));
   }
    public String comment() {
	return COMMENT;
    }
}
