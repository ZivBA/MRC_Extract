package meshi.applications.loopBuilding.applications;

import java.io.IOException;

import meshi.applications.prediction.GDTcalculator;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.MinimizerException;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.file.MeshiWriter;


/**
 * 
 * @author Nir
 *
 * This application decides whether to take the residue from the trivial homology or from 
 * the refined homology model (called here stupidly as 'ref' - but don't mix it with native). 
 * If the N,CA,C atoms of the trivial all fall within 'cutoff' of the corresponding atoms in the
 * refined, then the trivial is taken. Else, the refined residue is taken.
 * 
 * At the end, the model is superimposed on the NATIVE with a subset of residues within 1.5Ang. 
 * This is the only reason why we need the native structure here.
 *  
 */

public class MergeTempleteAndRefine extends MeshiProgram implements Residues, AtomTypes{ 

	private static String targetNumber = null;  
	private static double cutoffVal = -999;  

	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		mergeTempAndRef("PreModels/"+targetNumber+"_homo.pdb",
				"PreModels/"+targetNumber+"_ref.pdb",
				"Models/"+targetNumber+".pdb",
				"Natives/"+targetNumber+".pdb",
				cutoffVal);

	} // Of main



	public static void mergeTempAndRef(String tempFileName, String refFileName, String outputFileName, String nativeFileName, double cutoff) {	
		AtomList atomList = new AtomList();
		Protein templ = new Protein((new AtomList(tempFileName)).noOXTFilter().backbone(), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		Protein ref = new Protein((new AtomList(refFileName)).noOXTFilter().backbone(), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		boolean[] goodTemp = new boolean[ref.residues().size()];

		for (int res=0 ; res<ref.residues().size() ; res++) {
			goodTemp[res] = true;
			if (!ref.residues().residueAt(res).dummy() && (ref.residues().residueAt(res).ca()!=null)) {
				if ((templ.residue(ref.residues().residueAt(res).number)==null) ||
						(templ.residue(ref.residues().residueAt(res).number).ca()==null)) {
					goodTemp[res] = false;
				}
				else { 
					if (distRefTemp(ref.residues().residueAt(res).atoms(),templ.residue(ref.residues().residueAt(res).number).atoms(), cutoff)) {
						goodTemp[res] = true;
					}
					else {
						goodTemp[res] = false;
					}
				}
			}
		}

		
		for (int res=0 ; res<ref.residues().size() ; res++) {
			if (!ref.residues().residueAt(res).dummy() && (ref.residues().residueAt(res).ca()!=null)) {
				if ((templ.residue(ref.residues().residueAt(res).number)==null) ||
						(templ.residue(ref.residues().residueAt(res).number).ca()==null)) {
					atomList.add(ref.residues().residueAt(res).atoms().backbone()); // Remove the .backbone() filter to import all the atoms. 
					System.out.println("Took from Refined residue:" + ref.residues().residueAt(res).number + "   Because not appearing in template.");
				}
				else { 
					if (goodTemp[res] && (goodTemp[res-1] || goodTemp[res+1])) {
						atomList.add(templ.residue(ref.residues().residueAt(res).number).atoms().backbone());
					}
					else {
						atomList.add(ref.residues().residueAt(res).atoms().backbone()); // Remove the .backbone() filter to import all the atoms.
						System.out.println("Took from Refined residue:" + ref.residues().residueAt(res).number + "   Large dis.");
					}
				}
			}
		}
		Protein model = new Protein(atomList, new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		GDTcalculator.alignBySubset(new AtomList(nativeFileName), model.atoms(), 1.5);		

		// Outputting
		try {
			model.atoms().print(new MeshiWriter(outputFileName));
		} catch (IOException e) {
			System.out.println("\n\nWriting refined model to disk was unsuccessful\n\n");
		}
	}

	public static boolean distRefTemp(AtomList ref,AtomList templ, double cutoff) {
		int resNumber = ref.atomAt(0).residueNumber();
		if ((ref.findAtomInList("N", resNumber).distanceFrom(templ.findAtomInList("N", resNumber))<cutoff) &&
				(ref.findAtomInList("CA", resNumber).distanceFrom(templ.findAtomInList("CA", resNumber))<cutoff) &&
				(ref.findAtomInList("C", resNumber).distanceFrom(templ.findAtomInList("C", resNumber))<cutoff))  {
			return true;
		}
		else {
			return false;		
		}
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


		String errorMessage = ("\n                  ******************\n"+
				"Usage java -Xmx300m MergeTempleteAndRefine <Target number> <cutoff> \n"+
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));

		targetNumber = getOrderedArgument(args);
		if (targetNumber == null) throw new RuntimeException(errorMessage);
		System.out.println("# Target number is "+targetNumber);

		String cutoffString = getOrderedArgument(args);
		if (cutoffString == null) throw new RuntimeException(errorMessage);
		cutoffVal = Double.parseDouble(cutoffString);
		System.out.println("# cutoff is: "+cutoffVal);

		initRandom(999);
	}
}
