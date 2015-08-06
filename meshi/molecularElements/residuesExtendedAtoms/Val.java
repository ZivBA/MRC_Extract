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
 *       H   CB - CG1
 *           |
 *           CG2
 *
 **/
public class Val extends ResidueExtendedAtoms {
    public final Atom CG1, CG2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
 "                O\n"+
 "                |\n"+
 "       N - CA - C...n\n"+
 "       |   |\n"+
 "       H   CB - CG1\n"+
 "           |\n"+
 "           CG2\n";
    public Val(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public Val(Val old,int number){
         this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	 setSecondaryStructure(old.secondaryStructure());
	       setAccessibility(old.accessibility());
    }
    public Val(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(VAL, atomList, number, mode,  addAtomsFlag);

	Object[] temp = new Object[2];
	int i=0;
	temp[i++] = CG1 = getAtom("CG1",VCG1, atomList, this);
	temp[i++] = CG2 = getAtom("CG2",VCG2, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);

	if ((CB  != null)  && (CG1  != null)) bonds.add(CB.bond(CG1));
	if ((CB  != null)  && (CG2  != null)) bonds.add(CB.bond(CG2));
    }
    public String comment() {
	return COMMENT;
    }
}
