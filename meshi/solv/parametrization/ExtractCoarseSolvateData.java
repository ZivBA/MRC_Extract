package meshi.solv.parametrization;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.DecimalFormat;

import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.LennardJones.LennardJonesCreator;
import meshi.energy.compositeTorsions.ramachandran.RamachandranCreator;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.ResidueBuilder;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiAttribute;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;

/**
 * This program will collect statistics on the environments of coarse representaions: ROT1 and CBs of the MESHI atom types (190 of them). 
 * For each atom type we would wish to characterize its:
 * 1. Number of heavy atom neighbors when the model is in Rot1.
 * 2. Number of CBs in a certain cutoff.
 * 
 * @author Nir
 */
public class ExtractCoarseSolvateData extends MeshiProgram implements Residues, AtomTypes{ 

	// We get these user defined parameters from the command line
	// ----------------------------------------------------------
	// The structure database:
	static String listOfStructures = "C:/Users/Nir/Loop_Building_Project/listPISCES.txt";
	// Output file:
	static String outputFile = "All_Atom_Solvate_Data_5.3"; 
	// The carbon cutoff (carbon will counted in a sphere with this radius):
	static double cutoffNeighbor = 5.3;    	


	public static void main(String[] args){
		// ------------------------
		// User defined parameters:
		// ------------------------
		// The carbon bins:
//		double[] binsNeighbor = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200};
		double[] binsNeighbor = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592, 593, 594, 595, 596, 597, 598, 599, 600};
		// Minimum size protein:
		int min_Prot_Size = 79;    	
		// How many residues to ignore in locality
		int ignore_local = 1;
		// Maximal LJ energy per atom to consider
		double max_LJ_energy = 0.5;    	
		// CommandList:
		CommandList commands = new CommandList("commands");
		// ----------------------------------------

		// Local variables
		double dis;
		int ty,ty1,ty2;
		Atom atom,atom1,atom2;
		SolvateExtractionAttribute attribute,attribute1,attribute2;
		int indNeighbor;

		init(args); 

		// Creators for the energy terms, and separate energy terms
		EnergyCreator[] energyCreators = {
				new LennardJonesCreator(1.0)
		};
		EnergyCreator[] rottenCreators = {
				new RamachandranCreator(1.0)
		};
		
		TotalEnergy energy;

		// Dunbrack Lib
		DunbrackLib lib = new DunbrackLib(commands, 1.0, 2);

		// The data array:
		double[][] dataCBs;
		dataCBs = new double[Atom.numberOfTypes()][binsNeighbor.length];
		double[][] dataROT1;
		dataROT1 = new double[Atom.numberOfTypes()][binsNeighbor.length];

		// Going over the models
		String[] models = File2StringArray.f2a(listOfStructures);
		for (int i=0 ; i<models.length ; i++) { 		
			System.out.println("Reading: " + models[i]);

			// Creating the model
			Protein model = null;
			try {
			model = new Protein(new AtomList(models[i]), new ResidueExtendedAtoms(ADD_ATOMS));
			}
			catch (Exception e) {
				
			}
			if ((model != null) && (model.atoms().CAFilter().size()>min_Prot_Size)) {
				
				for (int res=0 ; res<model.residues().size(); res++) {
					if ((model.residues().residueAt(res).type<20) &&
							(model.residues().residueAt(res).type>-1))
						ResidueBuilder.buildCentroid(model.residues().residueAt(res));
				}

				
				// Excluding clashing residues				
				DistanceMatrix dm = new DistanceMatrix(model.atoms(), 5.5 , 1.0, 4);
				energy = new TotalEnergy(model, dm, rottenCreators, commands);
				energy.evaluateAtoms();
				double[] resRamachEnergies = new double[model.residues().size()];
				int validResidues = 0;
				for (int rr=0 ; rr<resRamachEnergies.length ; rr++) { 
					if (!model.residues().residueAt(rr).dummy()) {
						resRamachEnergies[rr] = model.residues().residueAt(rr).ca().energy();
						validResidues++;
					}
					else
						resRamachEnergies[rr] = -999999.0;
				}

				dm = new DistanceMatrix(model.atoms(), 5.5 , 1.0, 4);
				energy = new TotalEnergy(model, dm, energyCreators, commands);
				energy.evaluateAtoms();
				for (int c=0 ; c<model.atoms().size() ; c++) {
					atom = model.atoms().atomAt(c);
					if (atom.energy()>max_LJ_energy)
						for (int cc=0 ; cc<atom.residue().atoms().size(); cc++)
							atom.residue().atoms().atomAt(cc).addEnergy(100.0);
				}
				
				int rottenResidues = 0;
				for (int rr=0 ; rr<resRamachEnergies.length ; rr++) { 
					if (((model.residues().residueAt(rr).type==GLY) && (resRamachEnergies[rr]>1.65)) ||
						((model.residues().residueAt(rr).type!=GLY) && (resRamachEnergies[rr]>1.5))) {
						rottenResidues++;
						for (int c=0 ; c<model.residues().residueAt(rr).atoms().size() ; c++) {
							model.residues().residueAt(rr).atoms().atomAt(c).addEnergy(100.0);
						}						
					}
				}
				System.out.println("There were "+ rottenResidues + " rotten residues out of " + validResidues + "     (" + rottenResidues*1.0/validResidues + "%)");

				
				dm = new DistanceMatrix(model.atoms(), cutoffNeighbor+0.1 , 0.2, 4);
				for (int c=0; c<model.atoms().size() ; c++)
					model.atoms().atomAt(c).addAttribute(new SolvateExtractionAttribute());

				// ****************************************************
				// ************    CB CB CB CB CB  ********************
				// ****************************************************
				// Extracting the solvate data for the CB representaion
				DistanceList dl = dm.nonBondedList();
				for (int c=0 ; c<dl.size(); c++) {
					dis = dl.distanceAt(c).distance();
					if ((dis<cutoffNeighbor) && 
							(Math.abs(dl.distanceAt(c).atom1().residueNumber() - dl.distanceAt(c).atom2().residueNumber())>ignore_local)) {  
						atom1 = dl.distanceAt(c).atom1();
						atom2 = dl.distanceAt(c).atom2();
						ty1 = atom1.type;
						ty2 = atom2.type;
						attribute1 = (SolvateExtractionAttribute) atom1.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
						attribute2 = (SolvateExtractionAttribute) atom2.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
						if (atom2.name.equals("CB"))
							attribute1.addToSumOfNeighbors(1.0);
						if (atom1.name.equals("CB"))
							attribute2.addToSumOfNeighbors(1.0);
//						if (atom2.isBackbone && !atom2.isHydrogen)
//							attribute1.addToSumOfNeighbors(1.0);
//						if (atom1.isBackbone && !atom1.isHydrogen)
//							attribute2.addToSumOfNeighbors(1.0);
					}
				}

				// Add the protein's data to the entire database
				for (int c=0 ; c<model.atoms().size() ; c++) {
					atom = model.atoms().atomAt(c);
					ty = atom.type;
					attribute = (SolvateExtractionAttribute) atom.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
					if (atom.energy()<max_LJ_energy) {
						indNeighbor = findBinIndice(binsNeighbor, attribute.getSumOfNeighbors());
						dataCBs[ty][indNeighbor]++;
					}
					else {
						//System.out.println("Serious clash in: " + atom + "\n");
					}
				}


				// ****************************************************
				// ************    ROT1 ROT1 ROT1  ********************
				// ****************************************************
				// Extracting the solvate data for the ROT1
//				RotamericTools.putIntoRot1(model, dm, lib);
//				dm = new DistanceMatrix(model.atoms(), cutoffNeighbor+0.1 , 0.2, 4);
//				for (int c=0; c<model.atoms().size() ; c++)
//					((SolvateExtractionAttribute) model.atoms().atomAt(c).getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE)).reset();		    
//				dl = dm.nonBondedList();
//				for (int c=0 ; c<dl.size(); c++) {
//					dis = dl.distanceAt(c).distance();
//					if ((dis<cutoffNeighbor) && 
//							(Math.abs(dl.distanceAt(c).atom1().residueNumber() - dl.distanceAt(c).atom2().residueNumber())>ignore_local)) {
//						atom1 = dl.distanceAt(c).atom1();
//						atom2 = dl.distanceAt(c).atom2();
//						ty1 = atom1.type;
//						ty2 = atom2.type;
//						attribute1 = (SolvateExtractionAttribute) atom1.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
//						attribute2 = (SolvateExtractionAttribute) atom2.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
//						if (!atom2.isHydrogen)
//							attribute1.addToSumOfNeighbors(1.0);
//						if (!atom1.isHydrogen)
//							attribute2.addToSumOfNeighbors(1.0);
//					}
//				}
//
//				// Add the protein's data to the entire database
//				for (int c=0 ; c<model.atoms().size() ; c++) {
//					atom = model.atoms().atomAt(c);
//					ty = atom.type;
//					attribute = (SolvateExtractionAttribute) atom.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
//					if (atom.energy()<max_LJ_energy) {
//						indNeighbor = findBinIndice(binsNeighbor, attribute.getSumOfNeighbors());
//						dataROT1[ty][indNeighbor]++;
//						if ((atom.type==KNZ) && (indNeighbor==0))
//							System.out.println(atom);
//					}
//					else {
//						//System.out.println("Serious clash in: " + atom + "\n");
//					}
//				}
			}
		}


		// Outputting CBs
		try{
			DecimalFormat fmt = new DecimalFormat("0.##");    
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile+".CB.txt"));
			for (int atomType=0; atomType<Atom.numberOfTypes() ; atomType++) {
				//System.out.println("## type:" + Atom.type(atomType));
				for (int cnc=0; cnc<binsNeighbor.length ; cnc++) {
					//System.out.print(dataCBs[atomType][cnc] + " ");
					bw.write(dataCBs[atomType][cnc] + " ");
				}
				//System.out.println();
				bw.write("\n");
			}
			bw.close();
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}

		// Outputting ROT1s
//		try{
//			DecimalFormat fmt = new DecimalFormat("0.##");    
//			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile+".ROT1.txt"));
//			for (int atomType=0; atomType<Atom.numberOfTypes() ; atomType++) {
//				//System.out.println("## type:" + Atom.type(atomType));
//				for (int cnc=0; cnc<binsNeighbor.length ; cnc++) {
//					//System.out.print(dataCBs[atomType][cnc] + " ");
//					bw.write(dataROT1[atomType][cnc] + " ");
//				}
//				//System.out.println();
//				bw.write("\n");
//			}
//			bw.close();
//		}
//		catch(Exception e) {
//			throw new RuntimeException(e.getMessage());
//		}

	}


	public static int findBinIndice(double[] bins, double value) {
		double minDiff = Double.MAX_VALUE;
		int minInd = -1;
		for (int c=0; c<bins.length ; c++) {
			if (Math.abs(value-bins[c])<minDiff) {
				minDiff = Math.abs(value-bins[c]);
				minInd = c;
			}
		}
		return minInd;
	}


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
				"Usage java -Xmx1000m ExtractSolvateData <list of structures> <output file name> <the distance cutoff of CNC>\n"+
		"                    ******************\n");

		listOfStructures = getOrderedArgument(args);
		if (listOfStructures == null) throw new RuntimeException(errorMessage);
		System.out.println("# The structures are taken from: "+listOfStructures);

		outputFile = getOrderedArgument(args);
		if (outputFile == null) throw new RuntimeException(errorMessage);
		System.out.println("# Output will be written to: "+outputFile);

		String cutoffString = getOrderedArgument(args);
		if (cutoffString== null) throw new RuntimeException(errorMessage);
		cutoffNeighbor = (new Double(cutoffString)).doubleValue();
		System.out.println("# CNC cutoff value: "+ cutoffNeighbor);

		initRandom(333);
	}
}
