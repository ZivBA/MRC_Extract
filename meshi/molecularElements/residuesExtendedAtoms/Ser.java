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
 *          /
 *        OG
 **/
public class Ser extends ResidueExtendedAtoms {
    public final Atom OG;
    public static final String COMMENT = "From Levitt, JMB 168:592 (1983) table 2.\n"+
	"                O\n"+
	"                |\n"+
	"       N - CA - C...n\n"+
	"           |\n"+
	"           CB\n"+
	"          /\n"+
	"        OG\n";
    public Ser(int number, int mode, double x, double y, double z) {
	this(number, getList(x, y, z), mode, ADD_ATOMS);
    }
    public Ser(Ser old,int number){
         this(number,old.atoms(),old.getMode(),old.getAddAtomsFlag());
	 setSecondaryStructure(old.secondaryStructure());
	       setAccessibility(old.accessibility());
    }
    public Ser(int number, AtomList atomList, int mode, int addAtomsFlag) {
	super(SER, atomList, number, mode,  addAtomsFlag);

	Object[] temp = new Object[1];
	int i=0;
	temp[i++] = OG = getAtom("OG",SOG, atomList, this);
	for (i = 0; i <temp.length; i++)
	    if (temp[i] != null) atoms.add(temp[i]);

	if ((CB  != null)   && (OG  != null)) bonds.add(CB.bond(OG));
    }
    public String comment() {
	return COMMENT;
    }
}

