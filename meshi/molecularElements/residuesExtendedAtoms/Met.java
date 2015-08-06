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
 *           CG - SD - CE
 **/
public class Met extends ResidueExtendedAtoms {
    public final Atom CG, SD, CE;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                 O\n"+
	"                 |\n"+
	"        N - CA - C...n\n"+
	"        |   |\n"+
	"        H   CB\n"+
	"            |\n"+
	"            CG - SD - CE\n";
    public Met(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public Met(Met old,int number){
	this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	setSecondaryStructure(old.secondaryStructure());
	      setAccessibility(old.accessibility());
    }
    public Met(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(MET, atomList, number, mode,  addAtomsFlag);

	Object[] temp = new Object[3];
	int i=0;
	temp[i++] = CG = getAtom("CG",MCG, atomList, this);
	temp[i++] = SD = getAtom("SD",MSD, atomList, this);
	temp[i++] = CE = getAtom("CE",MCE, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);

	if ((CB  != null)   && (CG  != null)) bonds.add(CB.bond(CG));
	if ((CG  != null)   && (SD  != null)) bonds.add(CG.bond(SD));
	if ((SD  != null)   && (CE  != null)) bonds.add(SD.bond(CE));
    }
    public String comment() {
	return COMMENT;
    }
}
