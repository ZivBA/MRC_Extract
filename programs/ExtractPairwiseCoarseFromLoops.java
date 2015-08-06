package programs;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.DecimalFormat;

import meshi.geometry.DistanceMatrix;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.MinimizerException;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;
import meshi.util.rotamericTools.RotamericTools;

/**
 *<pre>
 * This program will minimize a batch of loops given in a list, and will compare the results to a reference.
 * Only the loop residues are compared in the RMS. There is no superposition, and it is assumed the rest of
 * the protein is aligned with the reference.
 * Unix usage:
 *     java -Xmx300m MinimizeBatchOfLoops <commands file name> <file with list of pdbs> <ref pdb file name> <loop starting resisue> <loop ending residue> <Wrg> <Wev> <Wsolv> <Whb> <Wprop> <Wramach> <output PDB extension string>
 *
 * <commands file name> - A text file containing the different flags and parameters required for 
 *                        the run.
 * <file with list of pdbs> - These will be minimized
 *
 * <ref pdb file name> - GDT and RMS will be calculated to this structure.
 * 
 * <W...> - The weights
 * 
 * <output PDB extension string> - Each minimized PDB will be saved with this string folowing its original name 
 *
 **/

public class ExtractPairwiseCoarseFromLoops extends MeshiProgram implements Residues, AtomTypes{ 

	private static CommandList commands; 
	private static String commandsFileName = null;
	private static String modelsFileName = null;  
	private static String refFileName = null;  
	private static int resStart = -999;
	private static int resEnd = -999;
	private static int toTakeRef = 20;
	private static int toTakeCorrect = 5;
	private static double natCO = 1.5; //Ang
	


	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		DunbrackLib lib = new DunbrackLib(commands, 1.0, 1);
		Protein reference = null;
		Protein model = null;
		DistanceMatrix distanceMatrix = null;
		DecimalFormat fmt2 = new DecimalFormat("0.###");

		// The data array:
		double[][][] dataCBcorrect = new double[20][21][25];
		double[][][] dataROT1correct = new double[20][21][25];
		double[][][] dataCBref = new double[20][21][25];
		double[][][] dataROT1ref = new double[20][21][25];
		int[] takenRot1correct = new int[resEnd+1];
		int[] takenRot1ref = new int[resEnd+1];
		boolean[] correctRes = new boolean[resEnd+1];


		// Loading the reference and the background model	
		reference = new Protein((new AtomList(refFileName)).noOXTFilter(), new ResidueExtendedAtoms(ADD_ATOMS));
		RotamericTools.putIntoRot1(reference, new DistanceMatrix(reference.atoms(), 5.5, 2.0, 3), lib);
		for (int res=0 ; res<reference.residues().size() ; res++) {
			if (!reference.residues().residueAt(res).dummy())
				RotamericTools.putOinSCcenter(reference.residues().residueAt(res));
		}
		model = new Protein((new AtomList(refFileName)).noOXTFilter(), new ResidueExtendedAtoms(ADD_ATOMS));
		for (int cc=0 ; cc<model.atoms().size() ; cc++)
			model.atoms().atomAt(cc).setChain(" ");
		RotamericTools.putIntoRot1(model, new DistanceMatrix(model.atoms(), 5.5, 2.0, 3), lib);
		for (int res=0 ; res<model.residues().size() ; res++) {
			if (!model.residues().residueAt(res).dummy())
				RotamericTools.putOinSCcenter(model.residues().residueAt(res));
		}
	
		
		// Looping on all the models
		String[] models = File2StringArray.f2a(modelsFileName);
		boolean notFull = true;
		for (int i=0 ; (i<models.length) && notFull ; i++) {
			try {	
				// Reading the model
				System.out.println("Now doing: " + models[i]);
				AtomList loop = new AtomList(models[i]);
				for (int c=0; c<loop.size() ; c++) {
					Atom atom = loop.atomAt(c);
					model.atoms().findAtomInList(atom.name(), atom.residueNumber()).setXYZ(atom.x(), atom.y(), atom.z());
				}
				RotamericTools.putIntoRot1(model, new DistanceMatrix(model.atoms(), 5.5, 2.0, 3), lib);
				for (int res=0 ; res<model.residues().size() ; res++) {
					if (!model.residues().residueAt(res).dummy())
						RotamericTools.putOinSCcenter(model.residues().residueAt(res));
				}
				findCorrectRes(reference,model,correctRes);
				
				// The extraction 
				for (int res=resStart ; res<=resEnd ; res++) {
					if ((correctRes[res] && (takenRot1correct[res]<toTakeCorrect)) || 
							(!correctRes[res] && (takenRot1ref[res]<toTakeRef))) {
						Residue residue = model.residue(res);
						for (int c1=0 ; c1<residue.atoms().size(); c1++) {
							for (int c2=0 ; c2<model.atoms().size(); c2++) {
								Atom atom1 = residue.atoms().atomAt(c1);
								Atom atom2 = model.atoms().atomAt(c2);
								if ((atom2.residueNumber()<resStart) || (atom2.residueNumber()>resEnd) ||
										(atom1.number()>atom2.number())){
									double dis = dis2atoms(atom1, atom2);
									if ((dis<23.0) && 
											(Math.abs(atom1.residueNumber() - atom2.residueNumber())>1)) {  
										int type1 = atom1.residue().type;
										int type2 = atom2.residue().type;
										if (type1>type2) {									
											int tmp = type1;
											type1 = type2;
											type2 = tmp;
										}
										if (correctRes[res]) {
											if (atom2.name.equals("CB") && atom1.name.equals("CB"))
												dataCBcorrect[type1][type2][(int) Math.round(dis)]++;
											if (atom2.name.equals("CA") && atom1.name.equals("CB"))
												dataCBcorrect[atom1.residue().type][20][(int) Math.round(dis)]++;
											if (atom2.name.equals("CB") && atom1.name.equals("CA"))
												dataCBcorrect[atom2.residue().type][20][(int) Math.round(dis)]++;
											if (atom2.name.equals("O") && atom1.name.equals("O"))
												dataROT1correct[type1][type2][(int) Math.round(dis)]++;
											if (atom2.name.equals("CA") && atom1.name.equals("O"))
												dataROT1correct[atom1.residue().type][20][(int) Math.round(dis)]++;
											if (atom2.name.equals("O") && atom1.name.equals("CA"))
												dataROT1correct[atom2.residue().type][20][(int) Math.round(dis)]++;
										}
										else {
											if (atom2.name.equals("CB") && atom1.name.equals("CB"))
												dataCBref[type1][type2][(int) Math.round(dis)]++;
											if (atom2.name.equals("CA") && atom1.name.equals("CB"))
												dataCBref[atom1.residue().type][20][(int) Math.round(dis)]++;
											if (atom2.name.equals("CB") && atom1.name.equals("CA"))
												dataCBref[atom2.residue().type][20][(int) Math.round(dis)]++;
											if (atom2.name.equals("O") && atom1.name.equals("O"))
												dataROT1ref[type1][type2][(int) Math.round(dis)]++;
											if (atom2.name.equals("CA") && atom1.name.equals("O"))
												dataROT1ref[atom1.residue().type][20][(int) Math.round(dis)]++;
											if (atom2.name.equals("O") && atom1.name.equals("CA"))
												dataROT1ref[atom2.residue().type][20][(int) Math.round(dis)]++;								
										}
									}
								}
							}
						}
						if (correctRes[res])
							takenRot1correct[res]++;
						else
							takenRot1ref[res]++;
					}
				}
								
				notFull=false;
				for (int res=resStart ; res<=resEnd ; res++) {
					if ((takenRot1correct[res]<toTakeCorrect) | (takenRot1ref[res]<toTakeRef))
						notFull=true;
				}
			}
			catch (Exception e) {
				System.out.println("Evaluation was not successful.\n");
				e.printStackTrace();
				System.out.println("\nContinueing to next model.\n");
			}
		}

		// Outputting CBs
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(modelsFileName+".nat.CB.txt"));
			for (int atomType1=0; atomType1<20 ; atomType1++) {
				for (int atomType2=atomType1; atomType2<21 ; atomType2++) {
					bw.write(atomType1 + " " + atomType2 + " ");
					for (int cnc=0; cnc<25 ; cnc++) {
						bw.write(dataCBcorrect[atomType1][atomType2][cnc] + " ");
					}
					bw.write("\n");
				}
			}
			bw.close();
			bw = new BufferedWriter(new FileWriter(modelsFileName+".ref.CB.txt"));
			for (int atomType1=0; atomType1<20 ; atomType1++) {
				for (int atomType2=atomType1; atomType2<21 ; atomType2++) {
					bw.write(atomType1 + " " + atomType2 + " ");
					for (int cnc=0; cnc<25 ; cnc++) {
						bw.write(dataCBref[atomType1][atomType2][cnc] + " ");
					}
					bw.write("\n");
				}
			}
			bw.close();			
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}

		// Outputting ROT1s
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(modelsFileName+".nat.ROT1.txt"));
			for (int atomType1=0; atomType1<20 ; atomType1++) {
				for (int atomType2=atomType1; atomType2<21 ; atomType2++) {
					bw.write(atomType1 + " " + atomType2 + " ");
					for (int cnc=0; cnc<25 ; cnc++) {
						bw.write(dataROT1correct[atomType1][atomType2][cnc] + " ");
					}
					bw.write("\n");
				}
			}
			bw.close();
			bw = new BufferedWriter(new FileWriter(modelsFileName+".ref.ROT1.txt"));
			for (int atomType1=0; atomType1<20 ; atomType1++) {
				for (int atomType2=atomType1; atomType2<21 ; atomType2++) {
					bw.write(atomType1 + " " + atomType2 + " ");
					for (int cnc=0; cnc<25 ; cnc++) {
						bw.write(dataROT1ref[atomType1][atomType2][cnc] + " ");
					}
					bw.write("\n");
				}
			}
			bw.close();			
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}
		
	} // Of main


	private static double dis2atoms(Atom atom1, Atom atom2) {
		return Math.sqrt((atom1.x() - atom2.x())*
				(atom1.x() - atom2.x()) +
				(atom1.y() - atom2.y())*
				(atom1.y() - atom2.y()) + 
				(atom1.z() - atom2.z())*
				(atom1.z() - atom2.z()));
	}

	private static void	findCorrectRes(Protein reference,Protein model,boolean[] correctRes) {
		for (int res=resStart; res<=resEnd ; res++) {
			Atom caRef = reference.residue(res).ca();
			Atom caModel = model.residue(res).ca();
			Atom cbRef = reference.residue(res).atoms().findAtomInList("CB", res);
			Atom cbModel = model.residue(res).atoms().findAtomInList("CB", res);
			if (cbRef==null) {
				if (dis2atoms(caRef, caModel)<natCO)
					correctRes[res] = true;
				else
					correctRes[res] = false;
			}
			else {
				if ((dis2atoms(caRef, caModel)<natCO) && (dis2atoms(cbRef, cbModel)<natCO))
					correctRes[res] = true;
				else
					correctRes[res] = false;				
			}
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


		String line;
		String errorMessage = ("\n                  ******************\n"+
				"Usage java -Xmx300m MinimizeBatchOfCoarseLoops <commands file name> <file with list of pdbs> <file with list of phipsi data> <ref pdb file name> " +
				"<loop starting resisue> <loop ending residue> <Wev> <Whb> <Wtorval> <Wtether> <output PDB extension string>\n"+
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));
		commandsFileName = getOrderedArgument(args);
		if (commandsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# commandsFileName = "+commandsFileName);

		commands = new CommandList(commandsFileName);

		modelsFileName = getOrderedArgument(args);
		if (modelsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# initial model file name is "+modelsFileName);

		refFileName = getOrderedArgument(args);
		if (refFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# reference file name is "+refFileName);

		initRandom(999);

		String tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		resStart = (new Integer(tmpString)).intValue();
		System.out.println("# Starting residue is " + resStart);

		tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		resEnd = (new Integer(tmpString)).intValue();
		System.out.println("# Ending residue is " + resEnd);

	}
}
