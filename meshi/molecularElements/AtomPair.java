package meshi.molecularElements;

/**
 **/	
public class AtomPair implements Comparable{
 
    private int numberOfRenumberings;

    private Atom atom1;
    private Atom atom2;
    public final Atom atom1() {return atom1;}
    public final Atom atom2() {return atom2;}

    private int smallNumber;
    private int largeNumber;
    public final int smallNumber() {return smallNumber;}
    public final int largeNumber() {return largeNumber;}

    private int atom1Number;
    private int atom2Number;
    public final int atom1Number() {return atom1Number;}
    public final int  atom2Number() {return atom2Number;}
    
    private int smallType;
    private int largeType;
    public final int smallType() {return smallType;}
    public final int largeType() {return largeType;}

    private boolean withinRmax;
    public final boolean withinRmax() {return withinRmax;}

    private boolean imortal = false;
    public final boolean imortal() { return imortal;}
    public void setImortal() { imortal = true;}

    public static final int ATOM_PAIR_CAPACITY = 2;
    private AtomList atoms;
    public final AtomList atoms() {return atoms;}
    
    public  AtomPair(Atom atom1, Atom atom2) {
	withinRmax = true;
	this.atom1 = atom1;
	this.atom2 = atom2;
	atoms = new AtomList(ATOM_PAIR_CAPACITY);
	atoms.fastAdd(atom1);
	atoms.fastAdd(atom2);
	numberOfRenumberings = 0;
	renumber();
    }

    protected void renumber() {
	numberOfRenumberings++;
	atom1Number = atom1.number();
	atom2Number = atom2.number();
	if (atom1.number() > atom2.number()) {
	    largeNumber = atom1.number();
	    smallNumber = atom2.number();
	}
	else {
	    largeNumber = atom2.number();
	    smallNumber = atom1.number();
	}
	if (atom1.type > atom2.type) {
	    largeType = atom1.type;
	    smallType= atom2.type;
	}
	else {
	    largeType = atom2.type;
	    smallType = atom1.type;
	}
    }

    public Atom sharedAtom(AtomPair other) {
	if (other.atom1 == atom1) return atom1;
	if (other.atom1 == atom2) return atom1;
	if (other.atom2 == atom1) return atom2;
	if (other.atom2 == atom2) return atom2;
	return null;
    }
    
   public void pairOutRmax() {withinRmax = false;}
   public void pairWithinRmax() {withinRmax = true;}

    public boolean equals(Object obj) {
	if (! (obj instanceof AtomPair) ) return false;
	AtomPair other = ( AtomPair) obj;
	return (((atom1 == other.atom1) && (atom2 == other.atom2)) ||
		((atom1 == other.atom2) && (atom2 == other.atom1)));
    }

    public int compareTo(Object obj) {
	AtomPair other = (AtomPair) obj;
// 	System.out.println("dddddddd \n"+
// 			   this+"\n"+other+"\n"+(largeNumber > other.largeNumber)+" "+
// 			   (largeNumber < other.largeNumber)+
// 			   " "+(smallNumber > other.smallNumber)+" "+(smallNumber < other.smallNumber));
	if (largeNumber > other.largeNumber) return 1;
	if (largeNumber < other.largeNumber) return -1;
	if (smallNumber > other.smallNumber) return 1;
        if (smallNumber < other.smallNumber) return -1;
	return 0;
    }
	
    public String laconic (String prompt){ 
	return prompt+" AtomPair "+atom1.number()+" "+atom2.number();
    }

    public String toString (){ 
	return " AtomPair (distance "+atom1.distanceFrom(atom2)+")\n"+
	    "\tatom1   :"+atom1+"\n"+
	    "\tatom2   :"+atom2;
    }
    public boolean frozen() {
	return (atom1.frozen() & atom2.frozen());
    }
}
