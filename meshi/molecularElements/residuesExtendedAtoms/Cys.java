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
 *           SG 
 **/
public class Cys extends ResidueExtendedAtoms {
    public final Atom SG;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"       |   |\n"+
	"       H   CB\n"+
	"           |\n"+
	"          SG \n";
    public Cys(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public Cys(Cys old,int number){
         this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	 setSecondaryStructure(old.secondaryStructure());
	       setAccessibility(old.accessibility());
    }
    public Cys(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(CYS, atomList, number, mode,  addAtomsFlag);

	Object[] temp = new Object[1];
	int i=0;
	temp[i++] = SG = getAtom("SG",CSG, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);

 	if ((CB != null)   && (SG != null))bonds.add(CB.bond(SG));
   }
    public String comment() {
	return COMMENT;
    }
}
