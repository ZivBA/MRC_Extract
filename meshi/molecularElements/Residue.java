package meshi.molecularElements;

import java.util.Arrays;
import java.util.Iterator;

import meshi.parameters.MeshiPotential;
import meshi.parameters.Residues;
import meshi.util.MeshiList;
import meshi.util.filters.Filter;
import meshi.util.string.StringList;

/**
 * A generic residue.
 **/
public class Residue implements ResidueCreator, Comparable, MeshiPotential, Residues{
    

    /**
     * Residue name.
     **/
    public final String name;

    /**
     * Residue number. Note that with the current set of constructors,
     * it is the responsibility of the creating class to provide a meaningful number.
     **/
    public final int number;

    /** 
     * which of the 20 possible it is.
     **/ 
    public final int type;

    /** 
     * Atoms of this residue;
     **/
    protected AtomList atoms;

    /**
     * Co-valiant bonds in the residue.
     **/
    protected AtomPairList bonds;

    protected String chain;

    /**
     * First atom of the residue (typicaly N).
     * The one conncted to the previous residue.
     **/
    protected Atom head;

    /**
     * last atom of the residue (typicaly C).
     * The one conncted to the next residue.
     **/
    protected Atom tail;
 
    protected Atom prevAtom;
    protected Atom nextAtom;
    /**
     * A list of Residue names3 (three letter code).
     **/
    private static StringList names3 = new StringList();

    /**
     * A list of Residue names1 (one letter code).
     **/
    private static StringList names1 = new StringList();

    /**
     * The list of residue names is done.
     **/
    private static boolean namesDone = false;

    /**
     * HELIX, SHEET ect.
     * Acceptable values may be found in meshi/parameters/MeshiPotential. 
     **/ 
    private String secondaryStructure = null; 
    private String accessibility = null; 

    protected boolean dummy;
    public final int mode;
    
    //----------------------------------------------------- constructors --------------------------------------------
    /**
     * Dummy residue.
     **/
    public Residue() {name = "UNK"; number = -1; type = -1;chain = " "; dummy = false; mode = NORMAL;}

    public Residue(int number, String oneLetter) {
	this.number = number;
	chain = " ";
	type = Residue.type(oneLetter);
	name = Residue.nameThreeLetters(type);
	dummy = false;
	mode = NORMAL;
    }


    /**
     * A key for searches in residue lists.
     **/
    public Residue(int number) {name = "UNK"; this.number = number; type = -1; chain = " "; dummy = false; mode = NORMAL;} 
 
   public Residue(String name, int type, int number) {
       this(name, type, number, new AtomList(), NORMAL);
   }

   public Residue(String name, int type, int number, int mode) {
	this(name, type, number, new AtomList(), mode);
    }
    public Residue(String name, int type, int number, 
		   AtomList atomList) {
	this(name, type, number,  atomList, NORMAL);
    }

    public Residue(String name, int type, int number, 
		   AtomList atomList, int mode) {
	this.name = name;
	this.number = number;
	this.type = type;
	atoms = atomList;
	bonds = new AtomPairList();
	chain = " ";
	prevAtom = nextAtom = head = tail = null;
	atoms.setResidue(this);
	this.mode = mode;
    }

    //------------------------------------------------- methods -----------------------------------------------
    /**
     * Returns a list of the residue atoms.
     **/
    public AtomList atoms() {return atoms;}

    /**
     * Returns a list of the residue bonds.
     **/
    public AtomPairList bonds() {return bonds;}

    public String toString() {
	String attributes = "";
	if (secondaryStructure != null) attributes += " SS: "+secondaryStructure;
	if (accessibility != null) attributes += " Access: "+ accessibility;
	return name+"_"+number+attributes;
    }	    

    public int getMode() {return mode;}

    public String test() {
	String out = "Testing "+name+" "+comment()+"\n";
	
	Atom atom;
	Iterator atomIter = atoms.iterator();
	out += "atoms:\n";
	while ((atom = (Atom) atomIter.next()) != null) 
	    out += atom.comment()+"\n";
	
 	AtomPair bond;
 	Iterator bondIter = bonds.iterator();
 	out += "bonds:\n";
 	while ((bond = (AtomPair) bondIter.next()) != null) 
 	    out += "\t"+bond.atom1().name+"\t"+bond.atom2().name+"\n";
	return out;
    }
    public Atom head() {return head;}
    public Atom tail() {return tail;}
    public Atom prevAtom() {return prevAtom;}
    public Atom nextAtom() {return nextAtom;}
    protected void setPrevAtom(Atom atom) {
// 	if (prevAtom != null) 
// 	    throw new RuntimeException("prevAtom already set");
	prevAtom = atom;
    }
    protected void setNextAtom(Atom atom) {
// 	if (nextAtom != null) 
// 	    throw new RuntimeException("nextAtom already set");
	nextAtom = atom;
    }
    public String comment() {
	return "Residue";
    }



	
    /**
     * residue number - 1.
     **/
    public int position() { return number - 1;}

    /**
     * create a residue based on an atoms list.
     * creates a dummyResidue if the list is null or empty.
     **/
    public Residue create(AtomList atoms, int residueNumber, int mode) {
	if ((atoms == null) || (atoms.size() == 0)) return new DummyResidue(residueNumber);
	Atom atom = atoms.atomAt(0);
        return new Residue(atom.residueName(), -1,residueNumber, atoms, mode);
    } 

    public Residue create(String name, int residueNumber, int mode, double x, double y, double z) {
	return new Residue(name, -1,residueNumber, null, mode);
    }
	 
    public int compareTo(Object obj) {
	Residue other = (Residue) obj;
	if (number < other.number) return -1;
	if (number > other.number) return 1;
	return 0;
    }

    /**
     * Fetches an atom by its name.
     **/
    public Atom get(String name) {
	Iterator atomsIter = atoms.iterator();
	Atom atom;
	while ((atom = (Atom) atomsIter.next()) != null)
	    if (atom.name.equals(name)) return atom;
	return null;
    }

    public final String secondaryStructure() {
	return secondaryStructure;
    }

    public void  setSecondaryStructure(String secondaryStructure) {
	this.secondaryStructure = secondaryStructure;
    }

    public final String accessibility() {
	return accessibility;
    }

    public void  setAccessibility(String accessibility) {
	this.accessibility = accessibility;
    }

    //------------------------------------------------- Residue types -----------------------------------------
    /**
     * A list of residue types.
     **/
    private static MeshiList types = new MeshiList(new IsType()); //Type Objects  stored by their intValue
    private static Object[] typesArray; //Type objects  stored by associated string.

    /**
     * The list of atom types is done.
     **/
    private static boolean typesDone = false;

    private static class Type implements Comparable{
	public int type;
	public String name3, name1;
	public int ca;
	public int cb;
	
	public Type (String name3, String name1, int ca, int type) {
	    this.type = type;
	    this.name3 = name3;
	    this.name1 = name1;
	    this.ca = ca;
	    if ((name3 == "GLY") | (name3 == "HEL") | (name3 == "UNK"))cb = -1;
	    else cb = Atom.type(name1+"CB");
	}
	
	public int compareTo(Object obj) {
	    return name1.compareTo(((Type) obj).name1);
	}
	
	public boolean equalse(Object obj) {
	    return compareTo(obj) == 0;
	}
    }
    
    public static class IsType implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof Type);
	}
    }

    /**
     * Adds a new residue name to the names list.
     **/
    public static int addName(String name3, String name1,int ca) {
	if (namesDone) throw new RuntimeException("Residue name added after the list is done");
	Type type = new Type(name3,name1,ca,types.size());
	if (types.contains(type)) throw new RuntimeException(name3+" already exists in the names3 list");
	types.add(type);
	return type.type;
    }
    /**
     * Adds a new residue name to the names list, and declare the list done.
     **/
    public static int addName(String name3, String name1, int ca, String done) {
	int t;
	if (! done.equals("DONE")) 
	    throw new RuntimeException("weird parameter to addName(String name3, "+
				       "String name1, String done): "+done);
	t = addName(name3, name1, ca);	
	namesDone = true;
	typesArray = types.toArray();
	Arrays.sort(typesArray);
	return t;
    }

    /**
     * A convertor from three letter code to one letter code.
     **/
    public static String three2one(String name3) {
	int i ;
	for (i = 0; i < typesArray.length; i++) {
	    if (((Type)typesArray[i]).name3.equals(name3)) return ((Type) typesArray[i]).name1;
	}
	return null;
    }

    /**
     * A convertor from one letter code to three letter code.
     **/
    public static String one2three(char c) {
	return one2three(""+c);
    }

    /**
     * A convertor from one letter code to three letter code.
     **/
    public static String one2three(String name1) {
	int i ;
	for (i = 0; i < typesArray.length; i++) {
	    if (((Type)typesArray[i]).name1.equals(name1)) return ((Type) typesArray[i]).name3;
	}
	return null;
    }
    /**
     * A convertor from CA type to one letter code.
     **/
    public static int fromCA(int ca) {
	int i ;
	for (i = 0; i < typesArray.length; i++) {
	    if (((Type)typesArray[i]).ca == ca) return ((Type) typesArray[i]).type;
	}
	return -1;
    }

    /**
     * A convertor from residue type to CA type.
     **/
    public static int caOf(int type) {
	int i ;
	for (i = 0; i < typesArray.length; i++) {
	    if (((Type)typesArray[i]).type == type) return ((Type) typesArray[i]).ca;
	}
	return -1;
    }

    public int caOf() {return caOf(type);}

    /**
     * A convertor from residue type to CB type.
     **/
    public static int cbOf(int type) {
	if (type == GLY) return -1;
	int i ;
	for (i = 0; i < typesArray.length; i++) {
	    if (((Type)typesArray[i]).type == type) return ((Type) typesArray[i]).cb;
	}
	return -1;
    }

    /**
     * A convertor from one letter code to CA type.
     **/
    public static int caOf(String oneLetter) {
	int i ;
	for (i = 0; i < typesArray.length; i++) {
	    if (((Type)typesArray[i]).name1.equals(oneLetter)) return ((Type) typesArray[i]).ca;
	}
	return -1;
    }

    /**
     * A convertor from one letter code to CB type.
     **/
    public static int cbOf(String oneLetter) {
	int i ;
	for (i = 0; i < typesArray.length; i++) {
	    if (((Type)typesArray[i]).name1.equals(oneLetter)) return ((Type) typesArray[i]).cb;
	}
	return -1;
    }

    /** 
     * Converts residue type to name.
     **/
    public static String nameThreeLetters(int i) {
	return ((Type) types.elementAt(i)).name3;
    }

    public String nameThreeLetters() {
	return ((Type) types.elementAt(type)).name3;
    }

    /** 
     * Converts residue type to name.
     **/
    public static String nameOneLetter(int i) {
	if (i == -1) return "-";
	return ((Type) types.elementAt(i)).name1;
    }

    public String nameOneLetter() {
	if (type == -1) return three2one(name);
	return ((Type) types.elementAt(type)).name1;
    }

    public static int type(String name) {
	if (name == null) throw new RuntimeException("name parameter is null");
	if (name.length() == 3) 
	    return type(three2one(name));
	else {
	    int pos = Arrays.binarySearch(typesArray,new Type("",name,-1,-1));
	    if (pos == -1) throw new RuntimeException("Cannot find type "+name);
	    return ((Type) typesArray[pos]).type;
	}
    }
    public Atom getAtom(String atomName){
          Iterator atomsIter = atoms.iterator();
	Atom atom;

	while ((atom = (Atom) atomsIter.next()) != null) {
	    if (atom.name().equals(atomName))
            return atom;
	}
	return null;
    }

    public Atom ca() {
	    return getAtom("CA");
    }

    public Atom carboxylC(){
        return getAtom("C");
    }

    public Atom amideN() {
    	return getAtom("N");
    }

    public boolean dummy() {return dummy;}

    // ------------------------ for homology modeling -----------------------------
    public void setPrevNextAtomsToNull(){
        prevAtom = null;
        nextAtom = null;

    }
    public void initiateAtoms(){
        setPrevNextAtomsToNull();
        Iterator atomsIter = atoms.iterator();
	    Atom atom;
	    while ((atom = (Atom) atomsIter.next()) != null) {
	        atom.setResidueNumber(number);
            atom.emptyBonded();
	    }
    }


    public void setSecondaryStructure(char ss) {
        if (ss == 'H') secondaryStructure = "HELIX";
        if (ss == 'E') secondaryStructure = "SHEET";
        if (ss == 'C') secondaryStructure = "COIL";
    }

    public void setAccessibility(char ac) {
        if (ac == 'A') secondaryStructure = "ACCESSIBLE";
        if (ac == 'B') secondaryStructure = "BURIED";
    }

  public void defrost() {atoms.defrost();}

  public void freeze() {atoms.freeze();}

    public void setChain(String chain) {
	for (Iterator iter = atoms.iterator(); iter.hasNext();) {
	    Atom atom = (Atom) iter.next();
	    atom.setChain(chain);
	}
    }
}
    
