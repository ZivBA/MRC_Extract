package meshi.applications.loopBuilding.applications;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.MinimizerException;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;
import meshi.util.file.MeshiWriter;

public class CompleteMissingResidues extends MeshiProgram implements Residues, AtomTypes {

	private static String alignmentFile_init = null;  
	
	private static String modelFileName = null;

	private static String outFileName = null;
	
	private static int completeFrom_init = -1;

	private static int completeTo_init = -1;

	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 

		// Reading the alignment file
		String[] alignmentStrings = File2StringArray.f2a(alignmentFile_init);
		int firstResInQuery = (new Integer(alignmentStrings[1])).intValue();
		System.out.println("Alignment in query starts at: " + firstResInQuery);
		String queryAlignment = alignmentStrings[3].trim().replace("-", "");
		System.out.println("SEQ: " + queryAlignment);

		
		Protein query = addResiduesToProtein(new ExtendedAtomsProtein(modelFileName,DO_NOT_ADD_ATOMS) , 
				firstResInQuery, queryAlignment, completeFrom_init, completeTo_init);				 
		
    	// Making the segmod file
		String querySegmod = query.getSequence();
		while (querySegmod.charAt(0)=='-') {
			querySegmod = querySegmod.substring(1);
		}
		System.out.println(querySegmod);
		System.out.println("Cut under line for SEGMOD alignment file.\n***********************************************************");
		String segmodString = querySegmod.length() + " " + "Template_" + alignmentFile_init + " XXX\n";
    	for (int qCounter=0; qCounter<querySegmod.length() ; qCounter++) {
    		segmodString += " " + querySegmod.charAt(qCounter);
    		if ((qCounter % 30) == 29) {
    			segmodString += "\n";
    		}
    	}
    	if ((querySegmod.length() % 30) != 0) {
    		segmodString += "\n";
    	}
    	segmodString += querySegmod.length() + " " + "SEGMOD_" + alignmentFile_init + " YYY\n";
    	for (int qCounter=0; qCounter<querySegmod.length() ; qCounter++) {
    		segmodString += " " + querySegmod.charAt(qCounter);
    		if ((qCounter % 30) == 29) {
    			segmodString += "\n";
    		}
    	}
    	if ((querySegmod.length() % 30) != 0) {
    		segmodString += "\n";
    	}
		System.out.println(segmodString);
		System.out.println("***********************************************************");
		
		
		// Writing to file
		try {
			query.atoms().noOXTFilter().print(new MeshiWriter(outFileName));
		}
		catch (Exception e) {
			System.out.print("\nThere was a problem writing the output:\n" + e + "\n\nContinueing...\n\n");
		}		
	}
	
	
	public static Protein addResiduesToProtein(Protein oldModel , int firstResInQuery, String queryAlignment, int completeFrom, int completeTo) {
		// Building the protein instance of the trivial homology
		AtomList largerList = new AtomList();
		int resNumCounterQ = firstResInQuery-1;
		Atom lastAdded=oldModel.residue(oldModel.firstResidue()).ca();
		for (int position=0; position<queryAlignment.length() ; position++) {
			resNumCounterQ++;
			if ((oldModel.residue(resNumCounterQ)!=null) &&
					!oldModel.residue(resNumCounterQ).dummy()) { // Coping atoms from existing model
				for (int atomCounter=0 ; atomCounter<oldModel.residue(resNumCounterQ).atoms().size() ; atomCounter++) { 
					largerList.add(new Atom(oldModel.residue(resNumCounterQ).atoms().atomAt(atomCounter)));		
					lastAdded = oldModel.residue(resNumCounterQ).atoms().atomAt(atomCounter);
				}
			}
			else {
				if ((resNumCounterQ>=completeFrom)&&(resNumCounterQ<=completeTo)) {
					System.out.println("Added residue: "+Residue.one2three(queryAlignment.charAt(position))+" " +
							resNumCounterQ);
					largerList.add(new Atom(lastAdded.x(),lastAdded.y(),lastAdded.z(),1.0,"CA",
							Residue.one2three(queryAlignment.charAt(position)),
							resNumCounterQ,-1));
				}					
			}
		}
		Protein query = new Protein(largerList, new ResidueExtendedAtoms(ADD_ATOMS));
		for (int cc=0 ; cc<query.atoms().size() ; cc++)
			query.atoms().atomAt(cc).setChain("A");
		
		return query;
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
		int zvl = ALA; // force the reading of "meshi.parameters.Residues"
		zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"


		String errorMessage = ("\n                  ******************\n"+
				"Usage java -Xmx600m CompleteMissingResidues <alignment file name> <model filename> <output filename> <complete from> <complete to> \n"+
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));

		alignmentFile_init = getOrderedArgument(args);
		if (alignmentFile_init == null) throw new RuntimeException(errorMessage);
		System.out.println("# Alignment file name is: "+alignmentFile_init);

		modelFileName = getOrderedArgument(args);
		if (modelFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Completeing: "+modelFileName);
		
		outFileName = getOrderedArgument(args);
		if (outFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Output to: "+outFileName);

		String tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		completeFrom_init = (new Integer(tmp.trim())).intValue();
		System.out.println("# Complete from: "+completeFrom_init);

		tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		completeTo_init = (new Integer(tmp.trim())).intValue();
		System.out.println("# Complete to: "+completeTo_init);

		initRandom(999);
	}	

} // Of AddMissingResidues
