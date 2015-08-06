package meshi.molecularElements.residuesExtendedAtoms;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
/**
 *<pre>
 *ARGININ
 *
 *              O
 *              |
 *     N - CA - C...n
 *     |   |
 *     H   CB                 
 *         |                  
 *         CG - CD - NE - CZ - NH1
 *                        |
 *                        NH2
  **/
public class Arg extends ResidueExtendedAtoms {
    public final Atom CG, CD, NE, HE, CZ, NH1, NH2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"              O\n"+
	"              |\n"+
	"     N - CA - C...n\n"+
	"     |   |\n"+
	"     H   CB\n"+                 
	"         |\n"+                  
	"         CG - CD - NE - CZ - NH1\n"+
	"                   |    |\n"+
	"                   HE   NH2\n";
    public Arg(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public Arg(Arg old,int number){
         this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	 setSecondaryStructure(old.secondaryStructure());
	       setAccessibility(old.accessibility());
    }
    public Arg(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(ARG, atomList, number, mode,  addAtomsFlag);
	Object[] temp = new Object[7];
	int i=0;
	temp[i++] = CG = getAtom("CG",RCG, atomList, this);
	temp[i++] = CD = getAtom("CD",RCD, atomList, this);
	temp[i++] = NE = getAtom("NE",RNE, atomList, this);
	if (NE != null) 
	    temp[i++] = HE = getAtom("HE",RHE, atomList, this);
	else HE = null;
	temp[i++] = CZ = getAtom("CZ",RCZ, atomList, this);
	temp[i++] = NH1 = getAtom("NH1",RNH, atomList, this);
	temp[i++] = NH2 = getAtom("NH2",RNH, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);
	
	if ((CB != null) && (CG != null)) bonds.add(CB.bond(CG));
	if ((CG != null) && (CD != null)) bonds.add(CG.bond(CD));
	if ((CD != null) && (NE != null)) bonds.add(CD.bond(NE));
	if ((NE != null) && (CZ != null)) bonds.add(NE.bond(CZ));
	if ((NE != null) && (HE != null)) bonds.add(NE.bond(HE));
	if ((CZ != null) && (NH1 != null)) bonds.add(CZ.bond(NH1));
	if ((CZ != null) && (NH2 != null)) bonds.add(CZ.bond(NH2));
    }
    public String comment() {
	return COMMENT;
    }
}
