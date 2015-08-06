package meshi.applications.loopBuilding.applications;

import java.io.File;
import java.io.IOException;
import java.util.StringTokenizer;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;
import meshi.util.file.MeshiWriter;

public class FindGoodSegmentsWithLoops extends MeshiProgram implements Residues, AtomTypes {


	private static String modelFileName = null;
	private static String compareToFileName = null;
	private static String loopDirPath = null;
	private static String gapsFileName = null;
	private static double goodCutoff = 1.5;
	private static final int NumPredToConsider = 1;

	
	public static void main(String[] args) {
		init(args); 
		
		Protein model = new ExtendedAtomsProtein(modelFileName,DO_NOT_ADD_ATOMS);
		boolean[] goodRes = new boolean[model.residues().size()+1];
		for (int c=0; c<goodRes.length ; c++)
			goodRes[c] = false;
		Protein compareTo = new ExtendedAtomsProtein(compareToFileName,DO_NOT_ADD_ATOMS);
		String[] gaps = File2StringArray.f2a(gapsFileName);
		for (int gapC=0 ; gapC<gaps.length ; gapC++) {
			StringTokenizer st = new StringTokenizer(gaps[gapC]);
			int gapFrom = (new Integer(st.nextToken())).intValue();
			int gapTo = (new Integer(st.nextToken())).intValue();
			for (int predC=0 ; predC<NumPredToConsider ; predC++) {
				if ((new File(loopDirPath+"/Loop_"+gapFrom+"_"+gapTo+"/phase2/"+predC+".pdb")).exists()) {
					boolean[] goodInLoop = new boolean[gapTo+1];
					AtomList gapAtoms = new AtomList(loopDirPath+"/Loop_"+gapFrom+"_"+gapTo+"/phase2/"+predC+".pdb");
					// Finding consistent atoms between loop and the 'compareTo' protein
					for (int resInLoop=gapFrom ; resInLoop<=gapTo ; resInLoop++) {
						if (compareTo.atoms().findAtomInList("CA", resInLoop)==null) {
							goodInLoop[resInLoop] = false;
						}
						else {
							goodInLoop[resInLoop] = false;
							if (compareTo.atoms().findAtomInList("CB", resInLoop)!=null) { // NOT a Glycin
								if ((compareTo.atoms().findAtomInList("CA", resInLoop).distanceFrom(gapAtoms.findAtomInList("CA", resInLoop))<goodCutoff) &&
										(compareTo.atoms().findAtomInList("CB", resInLoop).distanceFrom(gapAtoms.findAtomInList("CB", resInLoop))<goodCutoff)) {
									goodInLoop[resInLoop] = true;
								}
								else {
									goodInLoop[resInLoop] = false;
								}
							}
							else {
								if (compareTo.atoms().findAtomInList("CA", resInLoop).distanceFrom(gapAtoms.findAtomInList("CA", resInLoop))<goodCutoff) {
									goodInLoop[resInLoop] = true;
								}
								else {
									goodInLoop[resInLoop] = false;
								}								
							}
						}
					}
					
					// Analyzing data - is there a consistent pattern between the loops and the 'compareTo' Protein.
					for (int  resInLoop=gapFrom+4 ; resInLoop<=(gapTo-4) ; resInLoop++) {
						if (goodInLoop[resInLoop-1]&&goodInLoop[resInLoop]&&goodInLoop[resInLoop+1]) {
							goodRes[resInLoop-1] = true;
							goodRes[resInLoop] = true;
							goodRes[resInLoop+1] = true;
							System.out.println("The residues " + (resInLoop-1) + " " + (resInLoop) + " " + (resInLoop+1) + " were indicated by model " + predC + " from loop: " + gapFrom + "_" + gapTo);
						}
					}					
				}
			}
		}
		
		AtomList finalList = new AtomList();
		for (int c=0 ; c<compareTo.atoms().size() ; c++) {
			if (goodRes[compareTo.atoms().atomAt(c).residueNumber()])
				finalList.add(compareTo.atoms().atomAt(c));
		}
		
		try {
			finalList.print(new MeshiWriter(modelFileName+".common.pdb"));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
			
	}


	/** ================================= init =========================================
	 *
	 *A static function for parsing of the command line arguments and assigning the 
	 *variables commandsFileName, modelFileName and randomNumberSeed with the right inputs. Note that this
	 *static method is using parsing functions such as getOrderedArguments that are defined in MeshiProgram
	 *that MinimizeProtein inherits.
	 **/

	protected static void init(String[] args) {

		String errorMessage = ("\n                  ******************\n"+
				"Usage java -Xmx300m MinimizeBatchOfCoarseLoops <commands file name> <file with list of pdbs> <file with list of phipsi data> <ref pdb file name> " +
				"<loop starting resisue> <loop ending residue> <Wev> <Whb> <Wtorval> <Wtether> <output PDB extension string>\n"+
		"                    ******************\n");

		int zvl = ALA; // force the reading of "meshi.parameters.Residues"
		zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"

		initRandom(999);
		
		modelFileName = getOrderedArgument(args);
		if (modelFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# The model file name is: "+modelFileName);

		compareToFileName = getOrderedArgument(args);
		if (compareToFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Finding good segments in this model: "+compareToFileName);
		
		loopDirPath = getOrderedArgument(args);
		if (loopDirPath == null) throw new RuntimeException(errorMessage);
		System.out.println("# The loop dir path is: "+loopDirPath);

		gapsFileName = getOrderedArgument(args);
		if (gapsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# The gaps file name: "+gapsFileName);

		String tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		goodCutoff = (new Double(tmpString)).doubleValue();
		System.out.println("# Good overlap is below: " + goodCutoff);
	}	

} // Of AddMissingResidues
