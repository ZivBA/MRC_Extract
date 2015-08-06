package meshi.molecularElements;
import java.io.PrintWriter;
import java.util.Iterator;

import javax.xml.stream.events.StartDocument;

import meshi.PDB.PdbLine;
import meshi.PDB.PdbLineATOM;
import meshi.PDB.PdbReader;
import meshi.sequences.AtomAlignment;
import meshi.util.MeshiException;
import meshi.util.MeshiIterator;
import meshi.util.MeshiList;
import meshi.util.Rms;
import meshi.util.file.MeshiLineReader;
import meshi.util.filters.Filter;
/**
 * A list of Atoms
 **/
public class AtomList extends MeshiList {
    public MeshiLineReader sourceFile = null;


    /**
     * An empty atom list
     **/
    public AtomList() {
        this(new IsAtom());
    }

    /**
     * An empty atom list with an appointed capacity
     **/
    public AtomList(int capacity) {
	super(new IsAtom(),capacity);
    }

    public AtomList(Atom first) {
        this();
	add(first);
    }
    
    /**
     * An empty atom list
     **/
     protected AtomList(Filter filter) {
        super(filter);
    }
    /**
     * An atom list based on a PDB file
     * In a multiple model file (such as NMR's) only model 1 is read.
     * Currently whenever there are two alternative atom positions present,
     * only position A is taken.
     **/
    public AtomList(MeshiLineReader file, Filter lineFilter){
	this();
	PdbLine line;
	boolean contRead = true;
        Atom atom;
	
	sourceFile = file;
	if (!(file instanceof PdbReader)) 
	    throw new MeshiException("Cannot construct AtomList from file: "+file+"\n"+
			   "Currently only PDB formated files are used.");
	while (((line = ((PdbReader) file).readPdbLine()) != null) && contRead) {
	    if (line.isAModel() & (line.getModel() != 1))
		contRead = false;
	    if (lineFilter.accept(line)) {
		try {
		    atom = new Atom(line);
//		}
//		catch (Exception ex) {
//		    throw new RuntimeException("Failed to parse:\n"+line+"\n"+"in "+file+"\n"+ex);
//		}
              if (atom.alternateLocation().equals("") || atom.alternateLocation().trim().equals("A")) {
                     add(new Atom(line));
	      }
	      setComment(file.name());

                }
                catch (Exception ex) {
                     //System.out.println("Failed to parse:\n"+line+"\n"+"in "+file+"\n"+ex);
                }



	    }
	} 
	try {
	    file.close();
	}
	catch (Exception ex) { throw new RuntimeException("Cannot close "+file+"\n"+ex);}
    }
    /**
     * An atom list based on a PDB file
     **/
    public AtomList(PdbReader file) {
		this(file, new PdbLineATOM());
    }
    public AtomList(String dir, String fileName) {
		this(new PdbReader(dir, fileName));
    }
    public AtomList(String fileName, Filter lineFilter) {
		this(new PdbReader(fileName), lineFilter);
    }
    public AtomList(String fileName) {
		this(fileName, new PdbLineATOM());
    }
    /**
     * An atom list based on a residue list
     **/
    public  AtomList(ResidueList residues) {
	this();
	Iterator residueIter = residues.iterator();
	Iterator atomIter;
	Residue residue;
	Atom atom;	
	while((residue = (Residue) residueIter.next()) != null) {
	    atomIter = residue.atoms().iterator();
	    while((atom = (Atom) atomIter.next()) != null) 
		add(atom);
	}
    }

    /**
     * Extracts the list elements that pass the filter.
     **/
    public AtomList filter(Filter filter) {
	return (AtomList) filter(filter, new AtomList());
    }
	
    // fastAdd is equal to add(Object element)  in MeshiList, but checking of type is excluded
    protected boolean fastAdd(Atom element) {
	if (size < capacity) {
	    internalArray[size] = element;
	    size++;
	    modCount++;
	    return true;
	}
	else {
	    capacity *= 2;
	    Object[] newArray = new Atom [capacity];
	    for (int i = 0; i < size; i++)
		newArray[i] = internalArray[i];
	    internalArray = newArray;
	    return fastAdd(element) ;
	}
    }
    
    /**
     * RMS deviation between this list of atom and some other list.
     **/
    public Rms rms(AtomList otherList) {
	return new Rms(new AtomAlignment(this, otherList));
    }

    
    public double getRms(AtomList otherList) {
	return rms(otherList).getRms();
    }
    
    
    /**
     * Returns the atom in position <b>i</b> of the list.
     **/
    public Atom atomAt(int i) { return (Atom) elementAt(i);}
    public static class IsAtom implements Filter {
	public boolean accept(Object obj) {
		return (obj instanceof Atom); 
        }
    }

    public int whereIs(Atom atom) {
	for (int i = 0; i < size; i++) {
	    if (atomAt(i).equals(atom)) return i;
	}
	return -1;
    }

    /**
     * of gyration of the list atoms.
     * SQRT{ SIGMA_OVER_ALL_ATOMS([(x-X)^2+(y-Y)^2+(z-Z)^2]/N)}<br>
     * Where (x,y,z) are the atom coordinates and (X,Y,Z) are the 
     * coordinates of the center of mass.
     **/
    public double radius() {
              return  Math.sqrt(radius2() );
    }

    public double radius2() {
	double radius;
        MeshiIterator atomIter = meshiIterator(); 
	double dx, dy, dz, d2;
	double cmx, cmy, cmz; // center of mass x, y and z
        Atom atom;
	radius = 0;
	cmx = cmy = cmz = 0.0;
	while((atom = (Atom) atomIter.next()) != null) { 
	    cmx += atom.x();
	    cmy += atom.y();
	    cmz += atom.z();
	}
	cmx /= size();
	cmy /= size();
	cmz /= size();
	atomIter = meshiIterator();
	while((atom = (Atom) atomIter.next()) != null) { 
		dx = cmx - atom.x();
		dy = cmy - atom.y();
		dz = cmz - atom.z();
		d2 = dx*dx + dy*dy + dz*dz;
		radius += d2;
	}
	return radius = radius / size();
    }



    /**
     * Prints the atoms of the list along with their class.
     **/
    public void printClass() {
	Atom atom;
	Iterator atoms = iterator();
	while((atom = (Atom) atoms.next()) != null)
	    System.out.println(atom);
    }

    public void toFile(String header,String filename) {
        try{
            PrintWriter pw = new PrintWriter(filename);
            if(header!=null){
                pw.println(header);
            }
            Atom atom;
            Iterator atoms = iterator();
            while((atom = (Atom) atoms.next()) != null){
                pw.println(atom);
            }
            pw.close();
        }catch(Exception e){
            e.printStackTrace();
        }
    }

    /**
     * Sets the parameter to be the residue of this atom.
     **/
    public void setResidue(Residue residue) {
	Iterator atoms = iterator();
	Atom atom;
	while ((atom = (Atom) atoms.next()) != null)
	    atom.setResidue(residue);
    }
    public Atom[] toArrayOfAtoms() {
	Atom[] out = new Atom[size()];
	for (int i = 0; i < size(); i++)
	    out[i] = atomAt(i);
	return out;
    }
    
    public MeshiLineReader sourceFile() { return sourceFile;}
    
    public boolean frozenAtomsExist() {
	Iterator atoms = iterator();
	Atom atom;
	while((atom = (Atom) atoms.next()) != null)
	    if (atom.frozen()) return true;
	return false;
    }

    public void freeze() {
	Iterator atoms = iterator();
	Atom atom;
	while((atom = (Atom) atoms.next()) != null)
	    atom.freeze();
    }
    public void freeze(Filter filter) {
	Iterator atoms = iterator();
	Atom atom;
	while((atom = (Atom) atoms.next()) != null)
	    if (filter.accept(atom)) atom.freeze();
    }
    public void freezeAtomsWithCoordinates() {
	Filter filter = new HasCoordinates();
	Iterator atoms = iterator();
	Atom atom;
	while((atom = (Atom) atoms.next()) != null)
	    if (filter.accept(atom)) atom.freeze();
    }
    
    private class HasCoordinates implements Filter {
	public boolean accept(Object obj) {
	    Atom atom = (Atom) obj;
	    return atom.coordinates() != null;
	}
    }

    public void defrost() {
	Iterator atoms = iterator();
	Atom atom;
	while((atom = (Atom) atoms.next()) != null)
	    atom.defrost();
    }
    public AtomList fullOccupancyFilter() {
	Iterator atoms = iterator();
	Atom atom;
	AtomList out = new AtomList();
	while ((atom = (Atom) atoms.next()) != null) 
	    if (atom.occupancyNumber() > 0.99) out.add(atom);
	return out;
    }

    public AtomList noOXTFilter() {
	Iterator atoms = iterator();
	Atom atom;
	AtomList out = new AtomList();
	while ((atom = (Atom) atoms.next()) != null) 
	    if (atom.name().trim().compareTo("OXT") != 0) out.add(atom);
	return out;
    }

    public AtomList sidechain() {
    	return (AtomList) filter(new SidechainFilter(), new AtomList());
    }

    public static class SidechainFilter implements Filter {
    	public boolean accept(Object obj) {
    		Atom atom = (Atom) obj;
    		if (atom.name().equals("N")) return false; 
    		if (atom.name().equals("H")) return false; 
    		if (atom.name().equals("CA")) return false; 
    		if (atom.name().equals("CB")) return false; 
    		if (atom.name().equals("C")) return false; 
    		if (atom.name().equals("O")) return false; 
    		return true;
    	}
    }

    
    public AtomList backbone() {
	return (AtomList) filter(new BackboneFilter(), new AtomList());
    }

    public static class BackboneFilter implements Filter {
	public boolean accept(Object obj) {
	    Atom atom = (Atom) obj;
	    if (atom.name().equals("N")) return true; 
	    if (atom.name().equals("H")) return true; 
	    if (atom.name().equals("CA")) return true; 
	    if (atom.name().equals("CB")) return true; 
	    if (atom.name().equals("C")) return true; 
	    if (atom.name().equals("O")) return true; 
	    return false;
	}
    }
        
    public static class ClassCaCbFilter implements Filter {
    	public boolean accept(Object obj) {
    	    Atom atom = (Atom) obj;
    	    if (atom.name().equals("CA") || atom.name().equals("CB")) return true; 
    	    return false;
    	}
        }    
    
    public static class ClassCAFilter implements Filter {
	public boolean accept(Object obj) {
	    Atom atom = (Atom) obj;
	    if (atom.name().equals("CA")) return true; 
	    return false;
	}
    }

    public static class ClassNCaCFilter implements Filter {
    	public boolean accept(Object obj) {
    	    Atom atom = (Atom) obj;
    	    if (atom.name().equals("CA")) return true; 
    	    if (atom.name().equals("N")) return true; 
    	    if (atom.name().equals("C")) return true; 
    	    return false;
    	}
        }
    
    
    public static class NoAlternativeLocationFilter implements Filter {
        public boolean accept(Object obj) {
            Atom atom = (Atom) obj;
            if (atom.alternateLocation().equals("")) return true;
            return false;
        }
    }

    public boolean sortable() {return false;}

    public String toString() { return "AtomList with "+size()+" atoms";}

    public AtomList CAFilter() {
        Iterator atoms = iterator();
        Atom atom;
        AtomList out = new AtomList();
	out.sourceFile = sourceFile;
        while ((atom = (Atom) atoms.next()) != null)
            if (atom.name().trim().compareTo("CA") ==0) out.add(atom);
        return out;
    }

    public AtomList chainFilter(String chainID) {
        Iterator atoms = iterator();
        Atom atom;
        AtomList out = new AtomList();
        out.sourceFile = sourceFile;
        while ((atom = (Atom) atoms.next()) != null)
            if (atom.chain().trim().equals(chainID.trim())) out.add(atom);
        return out;
    }
    
    public void renumber() {
	for (int i = 0; i < size(); i++)
	    atomAt(i).setNumber(i);
    }
    
    public AtomList duplicate() {
	AtomList out = new AtomList();
	Iterator atoms = iterator();
	Atom atom;
	
	while ((atom = (Atom) atoms.next()) != null)
	    out.add(new Atom(atom));
	return out;
    }

    public AtomList mirror() {
	AtomList out = new AtomList();
	Iterator atoms = iterator();
	Atom atom;
	
	//public Atom(double x, double y, double z,
	//	   String name, Residue residue, int type){
	while ((atom = (Atom) atoms.next()) != null)
	    out.add(new Atom(atom.x(), atom.y(), -1*atom.z(), atom.name(), atom.residue(), Atom.type(atom.type())));
	return out;
    }
	
     public AtomList frozenAtoms() {
       return (AtomList) filter(new IsFrozenFilter(), new AtomList());
     }

     public AtomList defrostedAtoms() {
       return (AtomList) filter(new IsDefrostedFilter(), new AtomList());
     }

     public static class IsDefrostedFilter implements Filter {
       public boolean accept(Object obj) {
           Atom atom = (Atom) obj;
           if (!atom.frozen()) return true;
           return false;
       }
     }

     public static class IsFrozenFilter implements Filter {
       public boolean accept(Object obj) {
           Atom atom = (Atom) obj;
           if (atom.frozen()) return true;
           return false;
       }
     }

     public static class NonHydrogen implements Filter {
       public boolean accept(Object obj) {
           Atom atom = (Atom) obj;
           return !(atom.isHydrogen);
       }
     }

     public static class CbFilter implements Filter {
       public boolean accept(Object obj) {
           Atom atom = (Atom) obj;
           return (atom.name().equals("CB"));
       }
     }

     public static class ChainA implements Filter {
         public boolean accept(Object obj) {
             Atom atom = (Atom) obj;
             return (atom.chain().trim().equals("A"));
         }
       }

     public static class ChainB implements Filter {
         public boolean accept(Object obj) {
             Atom atom = (Atom) obj;
             return (atom.chain().trim().equals("B"));
         }
       }

     public static class KCB_Filter implements Filter {
         public boolean accept(Object obj) {
             Atom atom = (Atom) obj;
             return (atom.name().equals("CB") && atom.residueName().equals("LYS"));
         }
       }

     public static class KCA_Filter implements Filter {
         public boolean accept(Object obj) {
             Atom atom = (Atom) obj;
             return (atom.name().equals("CA") && atom.residueName().equals("LYS"));
         }
       }
     
    /**
     * Returns true if the parameter has the same length and the same composition in terms of atom names.
     * Current implementation, rather inefficient.
     **/
    public boolean sameAtoms(AtomList other) {
	if (size != other.size) return false;
	
	for (int iAtom = 0; iAtom < size; iAtom++) {
	    boolean notFound =true;
	    Atom atom = (Atom) elementAt(iAtom);
	    for (int jAtom = 0; (jAtom < size) & notFound; jAtom++) {
		Atom otherAtom = (Atom) elementAt(jAtom);
		if (atom.name().equals(otherAtom.name())) notFound = false;
	    }
	    if (notFound) return false;
	}
	return true;
    }
	    
    public AtomList lowestEnergyAtoms(double threshold) {
	
	if (size < 10) throw new RuntimeException("very weird");
	double sum = 0;
	double sum2 = 0;
	for (Iterator atomIter = iterator(); atomIter.hasNext();) {
	    Atom atom = (Atom) atomIter.next();
	    double e = atom.energy();
	    sum += e;
	    sum2 += e*e;
	}
	double avg = sum/size;
	double avg2 = sum2/size;
	double std = Math.sqrt(avg2-avg*avg);
	
	AtomList newList = new AtomList();
	for (Iterator atomIter = iterator(); atomIter.hasNext();) {
	    Atom atom = (Atom) atomIter.next();
	    if (atom.energy() < avg+std*threshold) newList.add(new Atom(atom));
	}
	if ((std < Math.abs(avg)) | (std < 10)) 
	    return newList;
	else return newList.lowestEnergyAtoms(threshold);	    
    }
    
    public void smoothEnergy(int window) {
	if (window % 2 == 0) throw new RuntimeException("window must be even");
	double[] energies = new double[size];
	for (int i = 0; i <size; i++)
	    energies[i] = ((Atom) elementAt(i)).energy();
	for (int i = window/2+1; i <size-window/2; i++) {
	    double avg = 0;
	    for (int j = 0; j < window; j++)
		avg += energies[i-window/2+j];
	    avg /= window;
	    ((Atom) elementAt(i)).resetEnergy();
	    ((Atom) elementAt(i)).addEnergy(avg);
	}
    }
    
    /** Gives the index of similar atom in name and residue number and residue type.
     * They need not be the same instance. If not found, returns -1.
     **/
    public int findIndexInList(Atom atom) {
    	for (int c=0 ; c<size() ; c++) 
    		if (atomAt(c).name().equals(atom.name()) &&
    		atomAt(c).residueName().equals(atom.residueName()) &&
    		(atomAt(c).residueNumber() == atom.residueNumber()))
    			return c;
    	return -1;
    }

    /** Gives the pointer to atom in list with similar name and residue number.
	* If not found, returns null.
     **/
    public Atom findAtomInList(String name , int resNumber) {
    	for (int c=0 ; c<size() ; c++) 
    		if (atomAt(c).name().equals(name) &&
    		(atomAt(c).residueNumber() == resNumber))
    			return atomAt(c);
    	return null;
    }
    
	/**
	 * Returns the atom with the specified name from the list.
	 **/
	public Atom getAtom(String atomName){
		for (int c=0 ; c<size() ; c++) 
			if (atomAt(c).name.equals(atomName)) 
				return atomAt(c);
		return null;
	}
    
    
    /** Gives the pointer to atom in list with similar name (i.e. a B or C or D atom) and residue number.
	* If not found, returns null.
     **/
    public Atom findAtomInList(char name  , int resNumber) {
    	for (int c=0 ; c<size() ; c++) 
    		if ((atomAt(c).name().length()>1) && (atomAt(c).name().charAt(1) == name) &&
    		(atomAt(c).residueNumber() == resNumber))
    			return atomAt(c);
    	return null;
    }

    /** Gives the atom index in list with similar name, residue number and chain ID.
	* If not found, returns -1.
     **/
    public int findAtomInList(String name , String chainID, int resNumber) {
    	for (int c=0 ; c<size() ; c++) 
    		if (atomAt(c).name().equals(name) &&  atomAt(c).chain().equals(chainID) &&
    		(atomAt(c).residueNumber() == resNumber))
    			return c;
    	return -1;
    }

    /** Like 'findAtomInList' above but returning an atom.
     **/
    public Atom findAtomInListReturningAtom(String name , String chainID, int resNumber) {
    	int c = findAtomInList(name , chainID, resNumber);
    	if (c==-1)
    		return null;
    	else
    		return atomAt(c);
    }
    
    /**
     * return a list of atoms that comprise the residue of atom 'statingAtomInex' and further down the list.
     */
    public AtomList getNextResidue(int startingAtomInex) {
    	AtomList returnList = new AtomList();
    	int resNumber = atomAt(startingAtomInex).residueNumber();
    	for (int c=startingAtomInex ; (c<size()) && (atomAt(c).residueNumber()==resNumber) ; c++) {
    		returnList.add(atomAt(c));
    	}
    	return returnList;
    }
    
    /**
     * Set all the chains of all atoms in the list to the parameter. 
     **/
    public void setChain(String newChain) {
    	for (int c=0 ; c<size() ; c++) {
    		atomAt(c).setChain(newChain);
    	}
    }

    public void moveCMtoOrigin() {
        MeshiIterator atomIter = meshiIterator(); 
		double cmx, cmy, cmz; // center of mass x, y and z
        Atom atom;
		cmx = cmy = cmz = 0.0;
		while((atom = (Atom) atomIter.next()) != null) { 
	    	cmx += atom.x();
	    	cmy += atom.y();
	    	cmz += atom.z();
		}
		cmx /= size();
		cmy /= size();
		cmz /= size();
		atomIter = meshiIterator();
		while((atom = (Atom) atomIter.next()) != null) { 
			atom.setX(atom.x()-cmx);
			atom.setY(atom.y()-cmy);
			atom.setZ(atom.z()-cmz);
		}
    }

    // This method return a new atom list which is the same as the current list (same Atom instances),
    // Except for residue resNum which is changed to a new residue type newType. Note that only the CA,
    // C,O,N, and H (in case the original residue was not proline) are present in the new type. The CB 
    // is not copied. 
    public AtomList mutateTo(int resNum , int newType) {
    	AtomList newList = new AtomList();
    	Atom atom;
    	for (int c=0 ; c<size() ; c++) {
    		atom = (Atom) elementAt(c);
    		if (atom.residueNumber()!=resNum)
    			newList.add(atom);
    		else if (atom.name.equals("CA") || atom.name.equals("C") || atom.name.equals("N") || 
    		         (atom.name.equals("H") && (newType!=12 /*Proline*/)) || atom.name.equals("O")) 
newList.add(new Atom(atom.x() , atom.y() , atom.z() , atom.name(), Residue.nameThreeLetters(newType), resNum, -1));    		         	
  		} 
    	return newList;
    }

    /**
     * Sorting the atoms in the list according to residue numbers
     */
    public AtomList sortByResidueNumber() {
        AtomList out = new AtomList();
        int minResNum=Integer.MAX_VALUE;
        int maxResNum=Integer.MIN_VALUE;
        for (int c=0 ; c<size() ; c++) {
        	if (atomAt(c).residueNumber()<minResNum)
        		minResNum = atomAt(c).residueNumber();
        	if (atomAt(c).residueNumber()>maxResNum)
        		maxResNum = atomAt(c).residueNumber();        	
        }
        for (int res=minResNum ; res<=maxResNum ; res++) {
        	for (int c=0 ; c<size() ; c++) {
        		if (atomAt(c).residueNumber()==res)
        			out.add(atomAt(c));
        	}
        }
        return out;
    }
    
    
    /**
     * Returns a string of the CAs residue names (in the order they appear in the list) 
     * in one-letter code. 
     */
    public String getSequence() {
    	String out ="";
    	for (int c=0 ; c<size() ; c++) {
    		if (atomAt(c).name().equals("CA")) {
    			out += Residue.three2one(atomAt(c).residueName());    			
    		}
    	}
    	return out;
    }
    
    /**
     * This method will add 1000 to the residue numbers of any new chain.
     */
    public void multiChain2meshi() {
    	int deltaRes = 1000;
    	String prevChain = "XXX";
    	for (int i=0 ; i<size() ; i++) {
    		if (!atomAt(i).chain.equals(prevChain) && !prevChain.equals("XXX")) {
    			deltaRes += 5000;
    		}
    		atomAt(i).setResidueNumber(atomAt(i).residueNumber() + deltaRes);
    		prevChain = atomAt(i).chain();
    	}
    }

    /**
     * This method does the opposite of the one above.
     */
    public void meshi2multiChain() {
    	int deltaRes = 1000;
    	String prevChain = "XXX";
    	for (int i=0 ; i<size() ; i++) {
    		if (!atomAt(i).chain.equals(" ")) {
    			if (!atomAt(i).chain.equals(prevChain) && !prevChain.equals("XXX")) {
    				deltaRes += 5000;
    			}
    			atomAt(i).setResidueNumber(atomAt(i).residueNumber() - deltaRes);
    			prevChain = atomAt(i).chain();
    		}
    	}
    }

    
}

