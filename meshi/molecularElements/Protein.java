package meshi.molecularElements;
import java.util.Iterator;

import meshi.PDB.OneOfTwenty;
import meshi.PDB.PdbLineATOM;
import meshi.PDB.PdbReader;
import meshi.geometry.AngleList;
import meshi.geometry.TorsionList;
import meshi.sequences.ResidueSequence;
import meshi.sequences.Sequence;
import meshi.sequences.SequenceAlignment;
import meshi.sequences.SequenceAlignmentCell;
import meshi.sequences.SequenceAlignmentColumn;
import meshi.sequences.SequenceList;
import meshi.util.KeyWords;
import meshi.util.MeshiProgram;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;
import meshi.util.string.StringList;
import meshi.util.string.StringParser;
/**
 * A protein chain.
 **/
public class Protein  extends MeshiProgram implements KeyWords{
    /**
     * The protein name.
     **/
    protected String name = "unkownProtein"; 
    
    /**
     * Often many models of the same protein are generated. 
     **/
    protected Integer modelNumber; 
    /**
     * A list of the protein's residues. 
     **/
    protected ResidueList residues; 

    /** 
     * A list of the protein's atoms.
     **/
    protected AtomList atoms; 

    /**
      * A list of the protein's bonds.
     **/
    protected AtomPairList bonds = null;
    protected AngleList angles = null;
    protected TorsionList torsions = null;

    protected int firstResidueIndex = -1;

    private boolean verbose = false;
    //----------------------------------------  constructors -------------------------------------------

    public Protein() {
    } 
    /**
     * Builds a protein from AA sequence and a secondary structure sequence.
     **/
    public Protein(String sequence, String secondaryStructure, String name, ResidueCreator creator) { 
    	residues = new ResidueList(sequence, creator);
	Iterator residuesIter = residues.iterator();
	Residue residue;
	int i = 0;
	char[] secondaryStructureArray = secondaryStructure.toCharArray();
	while ((residue = (Residue) residuesIter.next()) != null) {
	    char ss = secondaryStructureArray[i];
	    if (ss == 'H') residue.setSecondaryStructure("HELIX");
	    else if (ss == 'E') residue.setSecondaryStructure("SHEET");
	    else if (ss == 'C') residue.setSecondaryStructure("COIL");
	    else if (ss == 'A') residue.setSecondaryStructure("ALL");
	    else throw new RuntimeException("Unrecognized secondary structure "+ss);
	    i++;
	}
	this.name = name;
	atoms = new AtomList(residues);
	bonds = new AtomPairList(residues);

    }   

    //------------
   /**
     * Builds a protein from a PDB formatted file.
     * Uses generic atoms and residues. 
     **/
    public Protein(String fileName) { 
    	this(new PdbReader(fileName));
    }   

    public Protein(String fileName, ResidueCreator creator) {
        this(new AtomList(fileName),creator);
    }

    //------------

    /**
     * Builds a protein from a PDB formatted file.
     * Uses generic atoms and residues. 
     **/
    public Protein(MeshiLineReader file) { 
	this(file,  new PdbLineATOM());
    }


    //------------


    /**
     * Builds a protein from a PDB formatted file and a line filter.
     * 
     **/ 
    public Protein(MeshiLineReader file, Filter filter) {
	name = getProteinName(file);
	modelNumber = getModelNumber(file);
	atoms = new AtomList(file,filter);
	residues = new ResidueList(atoms, new Residue());
    }

    //------------


    /**
     * Builds a protein from a list of atoms and a ResidueCreator.
     * Allows Forcefield specific atoms and residued. 
     **/ 
     public Protein(AtomList atomList, ResidueCreator creator) {
	name = getProteinName(atomList.sourceFile());
	modelNumber = getModelNumber(atomList.sourceFile());
	residues = new ResidueList(atomList.filter(new OneOfTwenty()), 
				   creator);
	atoms = new AtomList(residues);
	bonds = new AtomPairList(residues);
// 	bonds.print();
    }
   
    //-------------------------------------------- Methods ----------------------------------------
    public void updateAtomList() {
	atoms =  new AtomList(residues);
    }

    public void updateBondsList() {
	bonds = new AtomPairList(residues);
    }
    public AtomPairList bonds() {return bonds;}

    
    /**
     * Extract the protein name from the filename.
     * expects name.pdb or name.modelNumber.pdb or pdbname.ent
     **/
    public static String getProteinName(MeshiLineReader file) {	
	if (file == null) return "unKnown";
	return getProteinName(file.path());
    }

    public static String getProteinName(String pathString) {
	StringList path = StringParser.breakPath(pathString);      
	String fileName = path.stringAt(path.size()-1);
	StringList temp = StringParser.breakFileName(fileName);
	String name;
	if (temp.size() <= 2) name = temp.stringAt(0);
        // was 	else  name = temp.stringAt(0)+"."+temp.stringAt(1);
        else{
            name = "";
            for(int i=0;i<temp.size()-1;i++)
                name += temp.stringAt(i)+(i<temp.size()-2?".":"");
        }
	if (name.startsWith("pdb"))
	    name = name.substring(3);
	return(name);
    }

    //------------


    /**
     * Extracts the model number from the filename.
     * Expects file name with the format name.modelNumber.pdb 
     **/
    private Integer getModelNumber(MeshiLineReader file) {
	if (file == null) return new Integer(0);
	StringList path = (StringParser.breakPath(file.path()));      
	String fileName = path.stringAt(path.size()-1);
	if ((StringParser.breakFileName(fileName)).size() >= 3) {
	    try {
		return new Integer((StringParser.breakFileName(fileName)).stringAt(1));
	    }
	    catch(Exception e) { return new Integer(0);} 
	}	
	else return new Integer(0);
    }

    /**
     * The protein name.
     **/
    public String name() { return name;}
    public int modelNumber() { 
	if (modelNumber == null) return 0;
	return modelNumber.intValue();
    }
 
   /**
     * A list of the protein's atoms.
     **/
    public AtomList atoms() {return atoms;}

  /**
     * A list of the protein's residues. 
     **/
    public ResidueList residues() {return residues;}


    /**
     * Residue number of the first residue in the residues list.
     * In most cases 1 but in many other cases some n > 1 typicaly because 
     * the N-terminus is not present in the PDB file
     **/
    public int firstResidue() {
	for (int i = 1; i <= residues.size(); i++) {
	    if (! residue(i).name.equals("UKN") &&  !residue(i).dummy()) return i;
	}
	return -1000;
    }

    /**
     * Residue number of the last residue in the residues list.
     **/
    public int lastResidue() {
    	return ((Residue) residues.last()).number;
    }
    
    
    /**
     * Returns the residue.
     **/
    public Residue residue(int residueNumber) { 
	return residues.residue(residueNumber);
    }

    /**
     * An iterator over the residues of the protein.
     **/
    public Iterator residueIterator() {
	return residues.iterator();
    }

    /**
     * Allow all atoms to move.
     **/
    public void defrost() {atoms.defrost();}

    public void freeze() {atoms.freeze();}
    public void freeze(Filter filter) {atoms.freeze(filter);}
    public void freezeAtomsWithCoordinates() {
	atoms.freezeAtomsWithCoordinates();
    }
    public String toString() {return name;}

    /**
     * Returns the specified atom.
     **/
    public Atom getAtom(String residueName, int residueNumber, String atomName){
	Iterator residuesIter = residues.iterator();
	Residue residue;
	Iterator atoms;
	Atom atom;
	while ((residue = (Residue) residuesIter.next()) != null) {
	    if ((residue.number == residueNumber) & residue.name.equals(residueName)) {
		atoms = residue.atoms.iterator();
		while ((atom = (Atom) atoms.next()) != null)
		    if (atom.name.equals(atomName)) {
			return atom;
		    }
	    }
	}
	return null;
    } 


    public void setResidues(ResidueList residueList) {
	residues = residueList;
	atoms = new AtomList(residues);
	atoms.renumber();
	bonds = new AtomPairList(residues); 
    }
    public int firstResidueIndex() {return firstResidueIndex;}

    public void allYouWantToKnow() {
	System.out.println("##################################################################################\n"+
			   "     Anything you ever wanted to know about "+this+" and never dared to ask\n"+
			   "##################################################################################\n");
 	System.out.println("========= residues ===========");
	residues().print();
	System.out.println("========= atoms  ===========");
	atoms().print();
	System.out.println("========= bonds ===========");
	bonds().print();
    }

		      
    /**
     *A method to set the secondary structure of a protein.
     *The input SS string currently support only C,E,H letters + A for the ALL type (= every SS is
     *posible). This method works only if the number of letter in the SS string equals to the 
     *number of non-dummy residues. The assignment of SS is than sequential in the residue numbers.
     **/
    public void setSS(String SS) {
    	
	Iterator residueIterator = residues().iterator();
	Residue res;
	int cc = -1; // residue number 0 is always dummy
	int resNum = -99999999;
    
	while ((res = (Residue) residueIterator.next()) != null) {
	    if (!(res instanceof DummyResidue)) {
		// Checking that the residues in the list are ordered by number
		if (resNum == -99999999)
		    resNum = res.number;
		else {
		    if (resNum >= res.number)
			throw new RuntimeException("\n\nThe residues in the residue list are not" +
						   " sorted by ascending residue number\n\n");
		    else
			resNum = res.number;
		}
		// Assigning the SS
		if (cc == SS.length())
		    throw new RuntimeException("\n\nThe secondary structure string provided is not long enough. \n\n");    			   
		if (SS.charAt(cc) == 'H') res.setSecondaryStructure("HELIX");
		else if (SS.charAt(cc) == 'E') res.setSecondaryStructure("SHEET");
		else if (SS.charAt(cc) == 'C') res.setSecondaryStructure("COIL");
		else if (SS.charAt(cc) == 'A') res.setSecondaryStructure("ALL");
		else throw new RuntimeException("\n\nUnrecognized secondary structure. Found: "+SS.charAt(cc)+"\n\n");    			   
	    }
	    cc++;
	}
	//preint the new residue list.
	residues.print(5," %-20s ");
    }	

    public void setSS(SequenceAlignment alignment) {    	
	Iterator residueIterator = residues().iterator();
	Residue residue;
        char prevSS = 'C';
        char ss;

	for (Iterator columns = alignment.iterator(); columns.hasNext();) {
		SequenceAlignmentColumn  column = (SequenceAlignmentColumn) columns.next();
		SequenceAlignmentCell cell = (SequenceAlignmentCell) column.cell(0);
		cell.addAttribute(column);
	}
	SequenceList tempList = new SequenceList(alignment);
	Sequence sequenceFromSS =(Sequence) tempList.elementAt(0);
	SequenceAlignment newAlignment = SequenceAlignment.identityAlignment(sequence(),sequenceFromSS);
	while ((residue = (Residue) residueIterator.next()) != null) {
	    if (!residue.dummy()) {
		SequenceAlignmentColumn alignmentColumn = (SequenceAlignmentColumn) newAlignment.getColumn(0,residue.number);
                if (alignmentColumn.cell(1).gap()) ss = prevSS;
		else {
		   SequenceAlignmentCell tempCell = (SequenceAlignmentCell) alignmentColumn.cell(1);
		   SequenceAlignmentColumn tempColumn = (SequenceAlignmentColumn) tempCell.getAttribute();
		   ss = tempColumn.getChar(1);
		   prevSS = ss;
		}

		if (ss == 'H') residue.setSecondaryStructure("HELIX");
		else if (ss == 'E') residue.setSecondaryStructure("SHEET");
		else if (ss == 'C') residue.setSecondaryStructure("COIL");
		else if (ss == 'A') residue.setSecondaryStructure("ALL");
		else throw new RuntimeException("\n\nUnrecognized secondary structure. Found: "+ss+"\n\n");    			   
	    }
	}
	//print the new residue list.
	if (verbose)  residues.print(5," %-20s ");
    }	


    public String getSequence(){
	String ans="";
	Iterator residuesIter = residueIterator();
	Residue res;
	residuesIter.next(); //to skip the dummy residue at position 0
	while ((res = (Residue)residuesIter.next())!=null) {
	    if (res instanceof DummyResidue) ans += "-";
	    else ans += res.nameOneLetter();
                                                }
	return ans;
    }
    

    public void printAtomsToFile(String fileName)throws Exception{
	
	MeshiWriter mw = new MeshiWriter(fileName);
	atoms.print(mw);
	mw.close();
    }

   public void setName(String name){this.name=name;}
    

    public Residue residueAt(int residueNumber) {
	return residues.residueAt(residueNumber);
    }
    
    public Sequence sequence() {
	String comment = name;
	String sequence = residues.toString();
	return new ResidueSequence(sequence, comment);
    }
    
    public void printCaspFormat(String fileName , 
    							String target,  /* with four digits T0302 T0277 etc. */
    							String author,  /* xxxxx-xxxxx-xxxxx */
    							int modelNum,  /* usually 1-5 */
    							String parent,  /* 1qwe_A or 1qwe  */
    							String method   /* Some description */	) throws Exception {
    	MeshiWriter mw = new MeshiWriter(fileName);
    	mw.println("PFRMAT TS ");
    	mw.println("TARGET " + target);
    	mw.println("AUTHOR " + author);
    	mw.println("METHOD " + method);
    	mw.println("MODEL " + modelNum);
    	mw.println("PARENT " + parent);
    	for (int c=0 ; c<atoms.size() ; c++) 
    		mw.println(atoms.atomAt(c)  + " ");
    	mw.println("TER ");
    	mw.println("END ");
    	mw.close();
    }
    
    public void defrostInRadiusAroundRes(Residue res , double radius) {
    	for (int c=0 ; c<residues.size() ; c++) 
    		if (!residues.residueAt(c).dummy())
    			if (residues.residueAt(c).ca().distanceFrom(res.ca())<radius)
    				residues.residueAt(c).atoms().defrost();
    }

    public static int[] getSeqOfProt(Protein prot, int fromRes, int toRes) {
    	int[] output = new int[toRes-fromRes+1];
    	for (int c=fromRes ; c<=toRes ; c++) {
    		if ((prot.residue(c) == null) || (prot.residue(c).ca() == null)) {
//    			throw new RuntimeException("\nCurrent setting asked for type of residue " + c + "which doesn't exist in the protein\n");
    			System.out.println("Warning: padding with alanines!!!!!!!");
    			output[c-fromRes] = 0;
    		}
    		else {
    			output[c-fromRes] = prot.residue(c).type;
    		}
    	}
    	return output;	
    }
    
}
