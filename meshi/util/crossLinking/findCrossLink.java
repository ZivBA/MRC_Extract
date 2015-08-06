package meshi.util.crossLinking;

import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.dssp.DSSP;


public class findCrossLink extends MeshiProgram implements Residues , ResidueMasses { 
   
    private static String complexFile = null;
 
    private static String dsspFile = null;

    private static String leftChain = null;  
 
    private static String rightChain = null;  

    public static void main(String[] args) {
	init(args); 

	AtomList complex = new AtomList(complexFile);
	DSSP dssp = new DSSP(dsspFile);
	
	leftChain = leftChain.trim();
	AtomList preprot1 = new AtomList();
	for (int cc=0 ; cc<complex.size() ; cc++)
		if (complex.atomAt(cc).chain().trim().equals(leftChain))
			preprot1.add(complex.atomAt(cc));
	Protein prot1 = new Protein(preprot1, new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));

	rightChain = rightChain.trim();
	AtomList preprot2 = new AtomList();
	for (int cc=0 ; cc<complex.size() ; cc++)
		if (complex.atomAt(cc).chain().trim().equals(rightChain))
			preprot2.add(complex.atomAt(cc));
	Protein prot2 = new Protein(preprot2, new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));

	
//	for (int cc=0 ; cc<left.atoms().size() ; cc++)
//		left.atoms().atomAt(cc).setChain(leftChain);
//	for (int cc=0 ; cc<right.atoms().size() ; cc++)
//		right.atoms().atomAt(cc).setChain(rightChain);
	
	crossLinkingTools.findCandidates(prot1, prot2, dssp, leftChain, rightChain, 5.0, 24.0, 0.1);
	
//	System.out.println("\n\n\n");
//	int[] pep1 = {2,12,0,4,8,7,3};
//	int[] pep2 = {7,13,11,9,9,8,7,16,11,9,14};
//	System.out.println("MASS:" + (MW_CL + MStools.massOfPeptide(pep1) + MStools.massOfPeptide(pep2) + 2*MW_OH + 2*MW_H) + "\n");
}



    /** ================================= init =========================================
     *
     *A static function for parsing of the command line arguments and assigning the 
     *variables commandsFileName, modelFileName and randomNumberSeed with the right inputs. Note that this
     *static method is using parsing functions such as getOrderedArguments that are defined in MeshiProgram
     *that MinimizeProtein inherits.
     **/
     
    protected static void init(String[] args) {
 
	/**** NOTE *** the next two lines. Because of a BUG in the Java VM, the 
	 * interfaces "Residues" and "AtomTypes" are not loaded automatically when MinimizeProtein initialize. 
	 * For this purpose these two lines are crucial wherever these two interfaces are implemented. The user might 
	 * rightfully feel that these two lines are "black magic" programming, but happily to our knowledge this is 
	 * the only bizarre phenomenon we are aware of in meshi.
	 **/
    	
    int zvl = KNZ;

	String errorMessage = ("\n                  ******************\n"+
			       "Usage java  findCrossLink <complex PDB file> <complex DSSP file> <left Chain> <right chain>\n"+
			       "                    ******************\n");
			      
	complexFile = getOrderedArgument(args);
	if (complexFile == null) throw new RuntimeException(errorMessage);
	System.out.println("# Complex file is: "+complexFile);

	dsspFile = getOrderedArgument(args);
	if (dsspFile == null) throw new RuntimeException(errorMessage);
	System.out.println("# DSSP of complex is: "+dsspFile);

	leftChain = getOrderedArgument(args);
	if (leftChain == null) throw new RuntimeException(errorMessage);
	System.out.println("# Left chain is: "+leftChain);

	rightChain = getOrderedArgument(args);
	if (rightChain == null) throw new RuntimeException(errorMessage);
	System.out.println("# Right chain is: "+ rightChain);
	
	initRandom(999);
    }
}
