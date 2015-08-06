package meshi.molecularElements;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.StringTokenizer;

import meshi.PDB.PdbLineATOM;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.util.filters.Filter;
/**
 * Proteins storage
 */
public class ComplexProtein extends Protein{
    protected Collection proteins;
    public ComplexProtein(){
        modelNumber = new Integer(0);
        proteins = new LinkedList();
    }
    public ComplexProtein(Protein p1){
        this();
        add(p1);
    }
    public ComplexProtein(Protein p1,Protein p2){
        this(p1);
        add(p2);
    }
    public ComplexProtein(Protein p1,Protein p2,Protein p3){
        this(p1,p2);
        add(p3);

    }
    /**
     * @@@ added by Ohad Givaty
     * @param pdbFileName String
     *
     */
    public  ComplexProtein(String pdbFileName){
	this(pdbFileName, new PdbLineATOM());
    }

    public  ComplexProtein(String pdbFileName, Filter lineFilter){
        this();
       AtomList fullList =  new AtomList(pdbFileName, lineFilter);
       Iterator fullListIter = fullList.iterator();
 
      if (!fullListIter.hasNext()) throw new RuntimeException("No atoms in atomList "+fullList.comment());
      Atom first = (Atom) fullListIter.next();
      if (!fullListIter.hasNext()) throw new RuntimeException("Not enough atoms in atomList "+fullList.comment());

      AtomList atomList = new AtomList();
      atomList.add(first);
      Atom atom = first;
      while (fullListIter.hasNext()) {
	  while (fullListIter.hasNext() & (atom.chain().equals(first.chain()))) {
	      atom = (Atom) fullListIter.next();
	      atomList.add(atom);
	  }
	  Protein prot = new Protein(atomList,new ResidueExtendedAtoms(1));
	  prot.setName(getName(fullList.sourceFile().path(),first.chain()));
	  first = atom;
	  this.add(prot);
      }
   }
    
    private String getName(String path, String chain) {
	StringTokenizer tokenizer = new StringTokenizer(path,"/\\");
	String fileName = "veryWeird_LooksLikeAbugIn_ComplexProtein_getName";
	while (tokenizer.hasMoreTokens()) fileName = tokenizer.nextToken();
	tokenizer = new StringTokenizer(fileName,".");
	String newName = tokenizer.nextToken();
	if (!chain.equals(" ")) newName +=chain;
	return newName;
    }

	   
    public Protein extractChain(String chain) {
	for (Iterator iter = proteins.iterator(); iter.hasNext();) {
	    Protein protein = (Protein) iter.next();
	    if (protein.atoms().atomAt(0).chain().equals(chain)) return protein;
	}
	return null;
    }

    public void add(Protein p){
        if(isEmpty()) {
            name = p.name;
	    residues = new ResidueList();
	    atoms = new AtomList();
	    bonds = new AtomPairList();
	}
        else
            name += "_"+name;
	residues.add(p.residues);
	atoms.add(p.atoms);
	bonds.add(p.bonds);
	proteins.add(p);
    }
    public boolean isEmpty(){return proteins.isEmpty();}
}
