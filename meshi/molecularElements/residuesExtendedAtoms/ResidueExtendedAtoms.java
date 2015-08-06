package meshi.molecularElements.residuesExtendedAtoms;
import java.util.Iterator;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.DummyResidue;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueCreator;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
public class ResidueExtendedAtoms extends Residue implements Residues, AtomTypes, ResidueCreator {
    private static final double RADIUS = 2;
    public final Atom H, N, CA, CB, C, O, OXT;    
    private static final int Ns[] =  {AN,  CN,  DN,  EN,  FN,  GN,  HN,  IN,  KN,  LN,  
				      MN,  NN,  PN,  QN,  RN,  SN,  TN,  VN,  WN,  YN};

    private static final int Hs[] =  {AH,  CH,  DH, EH,  FH,  GH,  HH,  IH,  KH,  LH,  
				      MH,  NH,  -1, QH,  RH,  SH,  TH,  VH,  WH,  YH};

    private static final int CAs[] = {ACA, CCA, DCA, ECA, FCA, GCA, HCA, ICA, KCA, LCA, 
				      MCA, NCA, PCA, QCA, RCA, SCA, TCA, VCA, WCA, YCA};

    private static final int CBs[] = {ACB, CCB, DCB, ECB, FCB,  -1, HCB, ICB, KCB, LCB, 
				      MCB, NCB, PCB, QCB, RCB, SCB, TCB, VCB, WCB, YCB};
    private static final int Cs[] =  {AC, CC, DC, EC, FC, GC, HC, IC, KC, LC, 
				      MC, NC, PC, QC, RC, SC, TC, VC, WC, YC};
    private static final int Os[] =  {AO, CO, DO, EO, FO, GO, HO, IO, KO, LO, 
				      MO, NO, PO, QO, RO, SO, TO, VO, WO, YO};
    private int addAtomsFlag = ADD_NOT_SET;
    
    /**
     * <pre>
     * Use this constructor to instantiate a creator object.
     * The addAtomsFlag parameter accept one of ADD_ATOMS, ADD_ATOMS_AND_FREEZE, 
     * DO_NOT_ADD_ATOMS or ADD_HYDROGENS_AND_FREEZE defined in meshi.parameters.Residues.
     * this flag affects the way the create methods behave. 
     * Use ADD_ATOMS for  Residue create(String name, int residueNumber, int mode, 
     *                                   double x, double y, double z).
     * Use ADD_ATOMS, ADD_ATOMS_AND_FREEZE, DO_NOT_ADD_ATOMS or ADD_HYDROGENS_AND_FREEZE for 
     * Residue create(AtomList atomList, int residueNumber, int mode) .
     *
     **/
    public ResidueExtendedAtoms(){
	this(ADD_ATOMS);
    }
    public ResidueExtendedAtoms(int addAtomsFlag){
	super();
	this.addAtomsFlag = addAtomsFlag;
	H = N = CA = CB = C = O = OXT = null;
    }


    /**
     * Do not use this constructor to instantiate a creator object.
     **/
    public ResidueExtendedAtoms(int type, AtomList atomList, int number, int mode, int addAtomsFlag) {
	super(nameThreeLetters(type), type, number,mode);
	this.addAtomsFlag = addAtomsFlag;
	String name =  nameThreeLetters(type);
	if ((mode == NORMAL) || (mode == CTER)) {
	    N = getAtom("N",Ns[type], atomList, this);
	    if ((N == null) | (type == PRO)) H = null;
	    else H = getAtom("H",Hs[type], atomList, this);
	}
	else  {
	    H = null;
	    N = getAtom("N",TRN,atomList, this);
	}

	CA = getAtom("CA",CAs[type], atomList, this);
	if ((mode == NORMAL) || (mode == NTER)) {
	    C = getAtom("C",Cs[type], atomList, this);
	    O = getAtom("O",Os[type], atomList, this);
	    OXT = null;
	}
	else {
	    C = getAtom("C",TRC, atomList, this);
	    O = getAtom("O",TRO, atomList, this);
	    OXT = getAtom("OXT",TRO, atomList, this);
	}

	if (type != GLY) CB = getAtom("CB",CBs[type], atomList, this);
	else CB = null;

	addAtomsAndBonds();
    }
	
    /**
     * Do not use this constructor to instantiate a creator object.
     **/
    public ResidueExtendedAtoms(int type, Atom PrevCa, Atom Ca, Atom NextCa, int number) {
	super(nameThreeLetters(type), type, number);
	String name =  nameThreeLetters(type);

	N = getAtom("N",Ns[type], PrevCa, Ca, 0.66, this);
	if (type == PRO) H = null;
	else H = new Atom(N.x(), N.y(), N.z(), 1.5, "H", this,Hs[type]);
	
	CA = new Atom(this, CAs[type], Ca, "CA");

	C = getAtom("C",Cs[type], Ca, NextCa, 0.33, this);
	O = new Atom(C.x(), C.y(), C.z(), 1.5, "O", this,Os[type]);
	OXT = null;

	if (type != GLY) CB = new Atom(CA.x(), CA.y(), CA.z(), 1.5, "CB", this,CBs[type]);
	else CB = null;

	addAtomsAndBonds();
    }

    private void addAtomsAndBonds() {
	if (H != null) atoms.add(H);
	if (N != null) atoms.add(N);
	if (CA != null) atoms.add(CA);
	if (CB != null) atoms.add(CB);
	if (C != null) atoms.add(C);
	if (O != null) atoms.add(O);
	if (OXT != null) atoms.add(OXT);

	if ((H != null) & (N != null))   bonds.add(H.bond(N));
	if ((CA != null) & (N != null))  bonds.add(N.bond(CA));
	if ((CA != null) & (C != null))  bonds.add(CA.bond(C));
	if ((O != null) & (C != null))   bonds.add(C.bond(O));
	if ((CB != null) & (CA != null)) bonds.add(CA.bond(CB));
	if ((OXT != null) & (C != null)) bonds.add(C.bond(OXT));
    }
	
    public Atom getAtom(String name, int type, Atom Ca1, Atom Ca2, double fraction, Residue residue) {
	if (type < 0) return null;
	double centerX = Ca1.x()+fraction*(Ca2.x() - Ca1.x());
	double centerY = Ca1.y()+fraction*(Ca2.y() - Ca1.y());
	double centerZ = Ca1.z()+fraction*(Ca2.z() - Ca1.z());

	return new Atom(centerX, centerY, centerZ, 1, name, residue, type);
    }


    public Atom getAtom(String name, int type, AtomList atomList, Residue residue) {
	if (type < 0) return null;
	double centerX = atomList.atomAt(0).x();
	double centerY = atomList.atomAt(0).y();
	double centerZ = atomList.atomAt(0).z();
	Iterator atomIter = atomList.iterator();
	Atom atom, newAtom;
	while ((atom = (Atom) atomIter.next()) != null) {
	    String atomName = nameConverter(atom.name);
	    if (atomName.equals(name)) {
		newAtom = new Atom(residue, type, atom, name);
		if ((addAtomsFlag == ADD_ATOMS_AND_FREEZE) ||
		    (addAtomsFlag == ADD_HYDROGENS_AND_FREEZE)) { 
		    newAtom.freeze();
		}
		return newAtom; 
	    }
	}
	if ((addAtomsFlag == ADD_ATOMS) ||
	    (addAtomsFlag == ADD_ATOMS_AND_FREEZE)||
	    ((addAtomsFlag == ADD_HYDROGENS_AND_FREEZE) && hydrogen(type))) {
	    newAtom = new Atom(centerX, centerY, centerZ, RADIUS,
			       name, residue, type);
	    return newAtom;
	}
	return null;
    }
		
    public static boolean hydrogen(int type) {
	if (type == AH) return true;
	if (type == CH) return true;
	if (type == DH) return true;
	if (type == EH) return true;
	if (type == FH) return true;
	if (type == GH) return true;
	if (type == HH) return true;
	if (type == IH) return true;
	if (type == KH) return true;
	if (type == LH) return true;
	if (type == MH) return true;
	if (type == NH) return true;
	if (type == QH) return true;
	if (type == RH) return true;
	if (type == SH) return true;
	if (type == TH) return true;
	if (type == VH) return true;
	if (type == WH) return true;
	if (type == YH) return true;
	if (type == HHD) return true;
	if (type == NHD1) return true;
	if (type == NHD2) return true;
	if (type == QHE1) return true;
	if (type == QHE2) return true;
	if (type == WHE) return true;
	if (type == HHE) return true;
	if (type == RHE) return true;
	return false;
    }
	
    public Atom head() {return  C;}
    public Atom tail() {return  N;}

    
    public Residue create(String name, int residueNumber, int mode, double x, double y, double z) {
	if (addAtomsFlag != ADD_ATOMS) throw new RuntimeException("This method may be called only "+ 
								  "if  addAtomsFlag == ADD_ATOMS");
	if (name.equals("A")) return new Ala(residueNumber, mode, x, y, z);
	if (name.equals("C")) return new Cys(residueNumber, mode, x, y, z);
	if (name.equals("D")) return new Asp(residueNumber, mode, x, y, z);
	if (name.equals("E")) return new Glu(residueNumber, mode, x, y, z);
	if (name.equals("F")) return new Phe(residueNumber, mode, x, y, z);

	if (name.equals("G")) return new Gly(residueNumber, mode, x, y, z);
	if (name.equals("H")) return new His(residueNumber, mode, x, y, z);
	if (name.equals("I")) return new Ile(residueNumber, mode, x, y, z);
	if (name.equals("K")) return new Lys(residueNumber, mode, x, y, z);
	if (name.equals("L")) return new Leu(residueNumber, mode, x, y, z);

	if (name.equals("M")) return new Met(residueNumber, mode, x, y, z);
	if (name.equals("N")) return new Asn(residueNumber, mode, x, y, z);
	if (name.equals("P")) return new Pro(residueNumber, mode, x, y, z);
	if (name.equals("Q")) return new Gln(residueNumber, mode, x, y, z);
	if (name.equals("R")) return new Arg(residueNumber, mode, x, y, z);

	if (name.equals("S")) return new Ser(residueNumber, mode, x, y, z);
	if (name.equals("T")) return new Thr(residueNumber, mode, x, y, z);
	if (name.equals("V")) return new Val(residueNumber, mode, x, y, z);
	if (name.equals("W")) return new Trp(residueNumber, mode, x, y, z);
	if (name.equals("Y")) return new Tyr(residueNumber, mode, x, y, z);
	return null;
    }

    public Residue create(AtomList atomList, int residueNumber, int mode) {
	if ((addAtomsFlag != ADD_ATOMS) &
	    (addAtomsFlag != ADD_ATOMS_AND_FREEZE) &
	    (addAtomsFlag != ADD_HYDROGENS_AND_FREEZE) &
	    (addAtomsFlag != DO_NOT_ADD_ATOMS))
	    throw new RuntimeException("This method may be called only "+ 
				       "if  addAtomsFlag == DO_NOT_ADD_ATOMS ADD_ATOMS, "+
				       "ADD_ATOMS_AND_FREEZE or ADD_HYDROGENS_AND_FREEZE");
	int type;
	String name = atomList.atomAt(0).residueName();
	if (name.equals(nameThreeLetters(type = ALA))) return new Ala(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = CYS))) return new Cys(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = ASP))) return new Asp(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = GLU))) return new Glu(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = PHE))) return new Phe(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = GLY))) return new Gly(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = HIS))) return new His(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = ILE))) return new Ile(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = LYS))) return new Lys(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = LEU))) return new Leu(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = MET))) return new Met(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = ASN))) return new Asn(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = PRO))) return new Pro(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = GLN))) return new Gln(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = ARG))) return new Arg(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = SER))) return new Ser(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = THR))) return new Thr(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = VAL))) return new Val(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = TRP))) return new Trp(residueNumber, atomList, mode, addAtomsFlag);
 	if (name.equals(nameThreeLetters(type = TYR))) return new Tyr(residueNumber, atomList, mode, addAtomsFlag);
	return null;
    }
    private static String nameConverter(String name) {
	if (name.equals("1HD2")) return "HD21";
	if (name.equals("2HD2")) return "HD22";
	if (name.equals("1HE2")) return "HE21";
	if (name.equals("2HE2")) return "HE22";
	return name;
    }
    
    protected static AtomList getList(double x, double y, double z) {
	AtomList out = new AtomList();
	out.add(new Atom(x, y, z, "dummy", new DummyResidue(-1), -1));
	return out;
    }
    public final Atom ca() { return CA;}
    public final Atom cb() { return CB;}
    public final Atom c() { return C;}
    public final Atom o() { return O;}
    public final Atom n() { return N;}
    public final Atom h() { return H;}

    public final Atom amideN() {return N;}

    // ***this is Guy's code***

    //a)
    /*      AA replacing methods        */
     
    public Arg toArg()  {
	Arg out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof Arg)) return (Arg)this;
	out=new Arg(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Gly toGly()  {
	Gly out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof Gly)) return (Gly)this;
	    out=new Gly(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Ala toAla() {
	Ala out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof Ala)) return (Ala)this;
	out=new Ala(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }
    public Val toVal() {
	Val out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof  Val)) return ( Val)this;
	out=new  Val(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Leu toLeu()  {
	Leu out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof  Leu)) return ( Leu)this;
	    out=new  Leu(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Ile toIle() {
	Ile out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof  Ile)) return ( Ile)this;
	    out=new  Ile(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Met toMet() {
	Met out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof Met )) return (Met )this;
	out=new Met(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Phe toPhe() {
	Phe out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof  Phe)) return (Phe)this;
	    out=new  Phe(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Trp toTrp() {
	Trp out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof  Trp)) return (Trp)this;
	out=new  Trp(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Pro toPro() {
	Pro out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof  Pro)) return (Pro)this;
	out=new  Pro(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Ser toSer() {
	Ser out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof  Ser)) return (Ser)this;
	out=new  Ser(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Thr toThr() {
	Thr out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof  Thr)) return (Thr)this;
	out=new  Thr(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Cys toCys() {
	Cys out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof  Cys)) return (Cys)this;
	    out=new  Cys(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Tyr toTyr() {
	Tyr out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof  Tyr)) return (Tyr)this;
	out=new  Tyr(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Asn toAsn() {
	Asn out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof  Asn)) return (Asn)this;
	out=new  Asn(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Gln toGln() {
	Gln out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof  Gln)) return (Gln)this;
	    out=new  Gln(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Asp toAsp() {
	Asp out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof  Asp)) return (Asp)this;
	    out=new  Asp(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Glu toGlu() {
	Glu out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof  Glu)) return (Glu)this;
	out=new  Glu(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public Lys toLys() {
	Lys out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof   Lys)) return ( Lys)this;
	out=new   Lys(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    public His toHis() {
	His out ;
	if ((N.type != TRN) && (C.type != TRC) && (this instanceof   His)) return ( His)this;
	out=new   His(this.number, NORMAL, CA.x(), CA.y(), CA.z());
	genericToAA(out);
	return out;
    }

    //b)
    protected void genericToAA(ResidueExtendedAtoms residue) {
	residue.setSecondaryStructure(secondaryStructure());
	residue.setAccessibility(accessibility());
	setBackboneAtoms(residue,this);
	Atom[] oldList = this.getSortedAtoms();
	Atom[] newList = residue.getSortedAtoms(); 
	int newPtr;
	newPtr=mapMatchingAtoms(oldList,newList);
	Atom temp;

	newList=listAddN(newList,residue);
	newPtr++;
	newList=listAddO(newList,residue);
	newPtr++;
	newList=listAddC(newList,residue);
	newPtr++;

	while (newPtr<newList.length) {

	    temp=getNextAtom(newPtr,newList,residue);
	    setNewCoordinates(newList[newPtr],temp.x(),temp.y(),temp.z());
	    newList[newPtr].setReliability(0.0);
	    newPtr++;

	}
	residue.CA.setReliability(BACKBONE_RELIABILITY);
	residue.N.setReliability(BACKBONE_RELIABILITY);
	residue.O.setReliability(BACKBONE_RELIABILITY);
	residue.C.setReliability(BACKBONE_RELIABILITY);
	return;
    
    } //end method genericToAA

    //c)
    private static void setBackboneAtoms(ResidueExtendedAtoms residue, ResidueExtendedAtoms in){

        setNewCoordinates(residue.CA,  in.CA.x(),  in.CA.y(),in.CA.z()) ;
        setNewCoordinates(residue.C,  in.C.x(),  in.C.y(),in.C.z())  ;
        setNewCoordinates(residue.O,  in.O.x(),  in.O.y(),in.O.z())  ;
        setNewCoordinates(residue.N,  in.N.x(),  in.N.y(),in.N.z())  ;

        if (in.H!=null && !(residue instanceof Pro)){
	    setNewCoordinates(residue.H,  in.H.x(),  in.H.y(),in.H.z())  ;
	}
            

    }
    //d)
    private Atom[] getSortedAtoms()
    {

	int atomsCounter=0;
	int copyItr=0;
	AtomList aList=this.atoms();
	if (aList==null) return null;

	for (int j=0 ; j<aList.size(); j++)
	    if  (((aList.atomAt(j)).name()).length()>1  &&  !((aList.atomAt(j)).name()).equals("OXT")  && !(((aList.atomAt(j)).name()).charAt(0)=='H')) atomsCounter++;


	Atom[] sortedList= new Atom[atomsCounter];
	for (int j=0 ; j<aList.size(); j++)
	    {//copy the pointers of the relevant atoms to a new list
		if  (((aList.atomAt(j)).name()).length()>1   &&  !((aList.atomAt(j)).name()).equals("OXT")    && !(((aList.atomAt(j)).name()).charAt(0)=='H'))
		    {
			sortedList[copyItr]=aList.atomAt(j);
			copyItr++;
		    }
	    }
	sortList(sortedList);
	return sortedList;

    }  //end method getSortedAtoms

    //e)
    private static int mapMatchingAtoms(Atom[] oldList,Atom[] newList){
	int newDist=1, oldDist=1;
    	int oldPtr=0,newPtr=0;

    	while (newPtr<newList.length && oldPtr<oldList.length )
	    {
	
		if ( (newPtr>0 && changedDistance(newList,newPtr)) && oldDist>=newDist)  newDist++;
		if (oldPtr>0 && changedDistance(oldList,oldPtr))  oldDist++;

		if (oldDist>=newDist )
		    {
        		newList[newPtr].setX(oldList[oldPtr].x());
        		newList[newPtr].setY(oldList[oldPtr].y());
        		newList[newPtr].setZ(oldList[oldPtr].z());
       		 	newPtr++;
		    }

		oldPtr++;
	    } //end while
    	
    	return newPtr;
			
    }
	
    //f)
    private static Atom[] listAddN(Atom[] newList, ResidueExtendedAtoms residue) {
			
        Atom temp[]= new Atom[newList.length+1];
        temp[0]=residue.N;
        
        for (int i=1; i<temp.length;i++)
	    temp[i]=newList[i-1];
	return temp;
    }
    private static Atom[] listAddO(Atom[] newList, ResidueExtendedAtoms residue) {

        Atom temp[]= new Atom[newList.length+1];
        temp[0]=residue.O;
		
        for (int i=1; i<temp.length;i++)
	    temp[i]=newList[i-1];
	return temp;
    }
    private static Atom[] listAddC(Atom[] newList, ResidueExtendedAtoms residue) {

        Atom temp[]= new Atom[newList.length+1];
        temp[0]=residue.C;
		
        for (int i=1; i<temp.length;i++)
	    temp[i]=newList[i-1];
	return temp;
    }
    //g)
    private static void setNewCoordinates(Atom atom, double X, double Y,double Z){
	
	atom.setX(X);
	atom.setY(Y);
	atom.setZ(Z);

    }
		
		
		
    private static Atom	getNextAtom(int newPtr, Atom[] newList,ResidueExtendedAtoms residue){
		
	double centerX,centerY,centerZ;       
	boolean condition3;		
	Atom chosenAtom,result=null;
	int modelAtom;
   
        if (newPtr>3 &&  !changedDistance(newList,newPtr) ) modelAtom=newPtr-2;

        else modelAtom=newPtr-1;
         
        centerX=newList[modelAtom].x() ;
        centerY=  newList[modelAtom].y() ;
        centerZ= newList[modelAtom].z() ;
        double minDistance=1.45,maxDistance=1.55;
        int itr=0;
	double DistFromAlpha=0;;
	for (int i=0;i<1000;i++){

	    double newDistFromAlpha;
						          
     	    do{

		chosenAtom=createNewAtom(centerX,centerY,centerZ,residue,maxDistance);//new Atom(centerX,centerY,centerZ,1.55,"locationGen",residue,1);
          
		itr++;
		if (itr==3000) { minDistance-=0.1; itr=0;}
		if (itr==1000) { maxDistance+=0.1;}
    		condition3=(minDist(chosenAtom.x(),chosenAtom.y(),chosenAtom.z(),newList,newPtr)< minDistance);

	    } while (condition3);
          
	    newDistFromAlpha=distance(chosenAtom.x(),chosenAtom.y(),chosenAtom.z(),newList[3].x(),newList[3].y(),newList[3].z());
	    if (newDistFromAlpha>DistFromAlpha) {DistFromAlpha=newDistFromAlpha;result=chosenAtom;}

	}		  
          

	return result;
		
    }
    //h)
    private final static double BACKBONE_RELIABILITY=200.0; 
    
    //i)
    private void sortList(Atom[] aList)   //performs insertion sort
    {
	if (aList==null) return;
	for (int i=1; i<aList.length; i++)
	    insert(aList,i);

    }

    //j)

    private static void insert(Atom[] aList, int i)
    {
	Atom value=aList[i];
	while (i>0 && compareAtoms(aList[i-1],value)>0)
	    {
		aList[i]=aList[i-1];
		i=i-1;
	    }

	aList[i]=value;
    }

    //k)

    private static boolean changedDistance(Atom[] list, int ptr){

	String firstAtomName= list[ptr].name();
	String secondAtomName = list[ptr-1].name();
	if (firstAtomName.charAt(1)!=secondAtomName.charAt(1) ) return true;
	return false;

    }

    //l)

		
    private static Atom createNewAtom(double centerX,double centerY,double centerZ,ResidueExtendedAtoms residue,double radius){

	/*      double atomsMaxDistance=1.55;

	double xDirection= ((atomCreationRand.nextInt()%1000)+1);
	double yDirection= ((atomCreationRand.nextInt()%1000)+1);
	double zDirection= ((atomCreationRand.nextInt()%1000)+1);
	double sum= Math.sqrt(Math.pow(xDirection,2)+Math.pow(yDirection,2)+Math.pow(zDirection,2));
	xDirection/=sum;
	yDirection/=sum;
	zDirection/=sum;
	xDirection*= atomsMaxDistance;
	yDirection*= atomsMaxDistance;
	zDirection*= atomsMaxDistance;                                             */
	return new Atom(centerX,centerY,centerZ,radius,"random",residue,1);
	//      return new Atom(centerX+xDirection,centerY+yDirection,centerZ+zDirection,"random",residue,1);


    }
    //m)
    public int getAddAtomsFlag() {return addAtomsFlag;}
    //n)

    private static double minDist(double x,double y,double z,Atom[] list,int place)
    {
	double min=10000;
	for (int i=0; i<place; i++)
	    {
		double dist=  distance(x,y,z,list[i].x(),list[i].y(),list[i].z());
		if (min>dist) min=dist;
	    }
	return min;
    }
    //o)
    private static double distance(double x1,double y1,double z1,double x2,double y2,double z2)
    {
        return Math.sqrt((x1-x2)*(x1-x2)+(y1-y1)*(y1-y2)+(z1-z2)*(z1-z2));
    }
    //p)
    private static int compareAtoms(Atom a1,Atom a2)
    {
	if ((a1==null) || (a2==null)) return 0;
	if ((a1.name()).equals(a2.name())) return 0;

	char c1,c2;
	char[] cArray={'A','B','G','D','E','Z','H'};
	char[] dArray={'1','2','3'};
	c1=(a1.name()).charAt(1);
	c2=(a2.name()).charAt(1);

	if (c1==c2)  //same distance from CA, compare positions( 1 or 2 )
	    {
		for   (int i=0; i<dArray.length ; i++)
		    {
			if ((a1.name()).charAt(2)==dArray[i]) return -1;
			if ((a2.name()).charAt(2)==dArray[i]) return 1;
		    }
	    }

	else
	    {
		for   (int i=0; i<cArray.length ; i++)
		    {
			if (c1==cArray[i]) return -1;
			if (c2==cArray[i]) return 1;

		    }
	    }



	return 0;
    }//end method compareAtoms

}
    
   
