package programs;

import meshi.applications.prediction.GDTcalculator;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;

/**
 *<pre>
 * This program will minimize a batch of proteins given in a list, and will compare the results to a reference.
 * 
 * Unix usage:
 *     java -Xmx300m MinimizeBatchOfProteins <commands file name> <file with list of pdbs> <ref pdb file name> <Wrg> <Wev> <Wsolv> <Whb> <Wprop> <Wramach>
 *
 * <commands file name> - A text file containing the different flags and parameters required for 
 *                        the run.
 * <file with list of pdbs> - These will be minimized
 *
 * <ref pdb file name> - GDT and RMS will be calculated to this structure.
 * 
 * <W...> - The weights
 *
 **/

public class EvaluateForKB01 extends MeshiProgram implements Residues, AtomTypes{ 

	private static String modelsFileName = null;  
	private static String refFileName = null;  
	private static String pathFirstEval = "";  
	private static String pathSecondEval = "";  


	public static void main(String[] args) {
		init(args); 
		Protein reference = null;
		Protein model = null;
		boolean fullMatchRefModel;

		// Loading the reference	
		reference = new Protein(refFileName, new ResidueExtendedAtoms(ADD_ATOMS));
//		RotamericTools.correctNomenclature(reference);


		// Looping on all the models
		String[] models = File2StringArray.f2a(modelsFileName);
		for (int i=0 ; i<models.length ; i++) { 		
			// Reading the model
			System.out.println("Trying to evaluate: " + pathFirstEval+ "/" + models[i]);
			try {
				model = new Protein(pathFirstEval+ "/" + models[i], new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
				for (int cc=0 ; cc<model.atoms().size() ; cc++)
					model.atoms().atomAt(cc).setChain("A");
				model.defrost();
//				RotamericTools.correctNomenclature(model);
				if (getMatchingAtoms(reference, model.atoms())==null)
					fullMatchRefModel = false;
				else
					fullMatchRefModel = true;


				//RMS - of original submission 
				System.out.println("999999 " + i + " 0000 " + models[i]);
				if (fullMatchRefModel)
					System.out.println("999999 " + i + " 1111 " + reference.atoms().CAFilter().getRms(getMatchingAtoms(reference,model.atoms()).CAFilter()) + 
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),0.5,1.0,2.0,4.0) +	
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),1.0,2.0,4.0,8.0) + " " +	
							reference.atoms().noOXTFilter().filter(new AtomList.NonHydrogen()).getRms(getMatchingAtoms(reference,model.atoms()).noOXTFilter().filter(new AtomList.NonHydrogen()))
							+ " " + getPercentageOfMatchingAtoms(reference.atoms(), model.atoms()));
				else
					System.out.println("999999 " + i + " 1111 " + (-1.0) + 
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),0.5,1.0,2.0,4.0) +	
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),1.0,2.0,4.0,8.0) + " " +	
							(-1.0)
							+ " " + getPercentageOfMatchingAtoms(reference.atoms(), model.atoms()));

				/* Writing dummy values for the energy */ 
				String tmpStr = "999999 " + i + " 2222 0 0 0 0 0 0 0 0 0 0 0 0 0";
				System.out.println(tmpStr);
				
//				System.out.println("Trying to evaluate: " + pathSecondEval + "/" + models[i]+".gbsa_min.pdb");
//				System.out.println("Trying to evaluate: " + pathSecondEval + "/" + models[i]);
//				System.out.println("Trying to evaluate: " + pathSecondEval + "/" + models[i]+"_KB-nH_bf.pdb.noRottenMinimialSCMOD.min.pdb");
//				System.out.println("Trying to evaluate: " + pathSecondEval + "/" + models[i]+".noRottenMinimialSCMOD.min.pdb");
				System.out.println("Trying to evaluate: " + pathSecondEval + "/" + models[i]+"_KB-nH_bf.pdb");
//				model = new Protein(pathSecondEval + "/" + models[i]+".gbsa_min.pdb", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
//				model = new Protein(pathSecondEval + "/" + models[i], new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
//				model = new Protein(pathSecondEval + "/" + models[i]+"_KB-nH_bf.pdb.noRottenMinimialSCMOD.min.pdb", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
//				model = new Protein(pathSecondEval + "/" + models[i]+".noRottenMinimialSCMOD.min.pdb", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
				model = new Protein(pathSecondEval + "/" + models[i]+"_KB-nH_bf.pdb", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
				for (int cc=0 ; cc<model.atoms().size() ; cc++)
					model.atoms().atomAt(cc).setChain("A");
				model.defrost();
//				RotamericTools.correctNomenclature(model);
				//RMS - Of altered submission
				if (fullMatchRefModel)
					System.out.println("999999 " + i + " 3333 " + reference.atoms().CAFilter().getRms(getMatchingAtoms(reference,model.atoms()).CAFilter()) + 
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),0.5,1.0,2.0,4.0) +	
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),1.0,2.0,4.0,8.0) + " " +	
							reference.atoms().noOXTFilter().filter(new AtomList.NonHydrogen()).getRms(getMatchingAtoms(reference,model.atoms()).noOXTFilter().filter(new AtomList.NonHydrogen()))
							+ " " + getPercentageOfMatchingAtoms(reference.atoms(), model.atoms()));
				else
					System.out.println("999999 " + i + " 3333 " + (-1.0) + 
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),0.5,1.0,2.0,4.0) +	
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),1.0,2.0,4.0,8.0) + " " +	
							(-1.0)
							+ " " + getPercentageOfMatchingAtoms(reference.atoms(), model.atoms()));
				/* Writing dummy values for the energy */ 
				tmpStr = "999999 " + i + " 4444 0 0 0 0 0 0 0 0 0 0 0 0 0";
				System.out.println(tmpStr);
			}
			catch (Exception e) {
				System.out.println("SKIPPING: Evaluation was not successful: " + models[i]);
			}
		}
	} // Of main




private static AtomList getMatchingAtoms(Protein refProt , AtomList protAtoms) {
	AtomList result = new AtomList();
	for (int c=0 ; c<refProt.residues().size() ; c++) 
		if (refProt.residues().residueAt(c).ca()!=null)
			if (protAtoms.findAtomInList("CA" , refProt.residues().residueAt(c).ca().residueNumber())==null) {
				System.out.println("MISMATCH: this atom is not found in model: " 
						+ refProt.residues().residueAt(c).ca());
				return null;
			}
	for (int c=0 ; c<protAtoms.size() ; c++)
		if (refProt.atoms().findAtomInList("CA" , protAtoms.atomAt(c).residueNumber())!=null)
			result.add(protAtoms.atomAt(c));

	return result;
}

private static double getPercentageOfMatchingAtoms(AtomList refList , AtomList modelList) {
	int n = 0;
	
	for (int c=0 ; c<modelList.size() ; c++) {
		if (modelList.atomAt(c).name().equals("CA"))
			if (refList.findAtomInList("CA", modelList.atomAt(c).residueNumber()) != null)
				n++;
	}

	return n*100.0/(refList.CAFilter().size());
}



/** ================================= init =========================================
 *
 *A static function for parsing of the command line arguments and assigning the 
 *variables commandsFileName, modelFileName and randomNumberSeed with the right inputs. Note that this
 *static method is using parsing functions such as getOrderedArguments that are defined in MeshiProgram
 *that MinimizeProtein inherits.
 **/

private static void init(String[] args) {

	/**** NOTE *** the next two lines. Because of a BUG in the Java VM, the 
	 * interfaces "Residues" and "AtomTypes" are not loaded automatically when MinimizeProtein initialize. 
	 * For this purpose these two lines are crucial wherever these two interfaces are implemented. The user might 
	 * rightfully feel that these two lines are "black magic" programming, but happily to our knowledge this is 
	 * the only bizarre phenomenon we are aware of in meshi.
	 **/
	int zvl = ALA; // force the reading of "meshi.parameters.Residues"
	zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"


	String line;
	String errorMessage = ("\n                  ******************\n"+
			"Usage java -Xmx300m MinimizeBatchOfProteins <commands file name> <file with list of pdbs> <ref pdb file name> " +
			"<Wrg>  <Wev> <Wsolv> <Whb> <Wprop> <Wramach> \n"+
	"                    ******************\n");

	if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));

	modelsFileName = getOrderedArgument(args);
	if (modelsFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# initial model file name is "+modelsFileName);

	refFileName = getOrderedArgument(args);
	if (refFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# reference file name is "+refFileName);

	initRandom(999);

	String tmpString = getOrderedArgument(args);
	if (tmpString!= null) { 
		pathFirstEval = new String(tmpString.trim());
		System.out.println("# pathFirstEval is " + pathFirstEval);
	}

	tmpString = getOrderedArgument(args);
	if (tmpString!= null) { 
		pathSecondEval = new String(tmpString.trim());
		System.out.println("# pathSecondEval is " + pathSecondEval);
	}
}
}
