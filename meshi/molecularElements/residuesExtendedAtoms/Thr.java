package meshi.molecularElements.residuesExtendedAtoms;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
/**
 *<pre>
 *                O
 *                |
 *       N - CA - C...n
 *           |
 *           CB
 *          /  \
 *        OG1  CG2
 **/
public class Thr extends ResidueExtendedAtoms{
    public final Atom OG1, CG2;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"           |\n"+
	"           CB\n"+
	"          /  \\\n"+
	"        OG1  CG2\n";
    public Thr(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public Thr(Thr old,int number){
         this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	 setSecondaryStructure(old.secondaryStructure());
	       setAccessibility(old.accessibility());
    }
    public Thr(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(THR, atomList, number, mode,  addAtomsFlag);

	Object[] temp = new Object[2];
	int i=0;
	temp[i++] = CG2 = getAtom("CG2",TCG, atomList, this);
	temp[i++] = OG1 = getAtom("OG1",TOG, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);

 	if ((CB  != null)   && (OG1 != null)) bonds.add(CB.bond(OG1));
 	if ((CB  != null)   && (CG2  != null)) bonds.add(CB.bond(CG2));
    }
    public String comment() {
	return COMMENT;
    }
}
