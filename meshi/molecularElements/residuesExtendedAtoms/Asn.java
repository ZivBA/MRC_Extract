package meshi.molecularElements.residuesExtendedAtoms;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
/**
 *<pre>
 *                 O
 *                |
 *       N - CA - C...n
 *       |   |
 *       H   CB
 *           |
 *           CG - ND2 - HD21
 *           |    |
 *           OD1  HD22
 **/
public class Asn extends ResidueExtendedAtoms {
    public final Atom CG, OD1, ND2, HD21, HD22;
    public static final String NAME = "ASN";
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                 O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       H   CB\n"+
	"           |\n"+
	"          CG - ND2 - HD21\n"+
	"           |    |\n"+
	"           OD1  HD22\n";
    public Asn(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public Asn(Asn old,int number){
         this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	 setSecondaryStructure(old.secondaryStructure());
	       setAccessibility(old.accessibility());
    }
   public Asn(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(ASN, atomList, number, mode,  addAtomsFlag);
	Object[] temp = new Object[5];
	int i=0;
	temp[i++] = CG = getAtom("CG",NCG, atomList, this);
	temp[i++] = ND2 = getAtom("ND2",NND, atomList, this);
	if (ND2 != null) {
	    temp[i++] = HD21 = getAtom("HD21",NHD1, atomList, this);
	    temp[i++] = HD22 = getAtom("HD22",NHD2, atomList, this);
	}
	else HD21 = HD22 = null;
	temp[i++] = OD1 = getAtom("OD1",NOD, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);
	
	if ((CB != null)   && (CG != null))  bonds.add(CB.bond(CG));
	if ((CG != null)   && (OD1 != null)) bonds.add(CG.bond(OD1));
	if ((CG != null)   && (ND2 != null)) bonds.add(CG.bond(ND2));
	if ((HD21 != null) && (ND2 != null)) bonds.add(ND2.bond(HD21));
	if ((HD22 != null) && (ND2 != null)) bonds.add(ND2.bond(HD22));
    }
    public String comment() {
	return COMMENT;
    }
}
