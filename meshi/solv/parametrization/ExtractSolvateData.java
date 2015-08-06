package meshi.solv.parametrization;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.DecimalFormat;

import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.LennardJones.LennardJonesCreator;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.hydrogenBond.AbstractHydrogenBond;
import meshi.geometry.hydrogenBond.AbstractHydrogenBondList;
import meshi.geometry.hydrogenBond.Dahiyat.DahiyatHighAccuracyParamaters;
import meshi.geometry.hydrogenBond.Dahiyat.DahiyatHydrogenBondListNoDuplications;
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
 * This program will collect statistics on the carbon environments and hydrogen bonding of the MESHI atom types (190 of them). 
 * For each atom type we would wish to characterize (at the same time) its:
 * 1. Number of heavy atom neighbors.
 * 2. Number of hydrogen bonds.
 * 3. Number of salt bridges.
 * 
 * @author Nir
 */
public class ExtractSolvateData extends MeshiProgram implements Residues, AtomTypes{ 

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
		double[] binsNeighbor = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200};
		// The hb bins:
		double[] binsHB = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5};
		// The salt bridge bins:
		double[] binsSB = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5};
		// Minimum size protein:
		int min_Prot_Size = 79;    	
		// Maximal LJ energy per atom to consider
		double max_LJ_energy = 0.5;    	
		// CommandList:
		CommandList commands = new CommandList("commands");
		// ----------------------------------------

		// Local variables
		double dis;
		Atom atom,atom1,atom2;
		SolvateExtractionAttribute attribute,attribute1,attribute2;
		AbstractHydrogenBond hb;
		int indNeighbor,indHB,indSB;

		init(args); 

		// Creators for the energy terms, and separate energy terms
		EnergyCreator[] energyCreators = {
				new LennardJonesCreator(1.0)
		};
		TotalEnergy energy;
		AbstractHydrogenBondList hbList;

		// The data array:
		double[][][][] data;
		data = new double[Atom.numberOfTypes()][binsNeighbor.length][binsHB.length][binsSB.length];

		// Going over the models
		String[] models = File2StringArray.f2a(listOfStructures);
		for (int i=0 ; i<models.length ; i++) { 		
			System.out.println("Reading: " + models[i]);

			/* If you need to add hydrogens to the models you will need the following commented code.
			 * --------------------------------------------------------------------------------------
			Protein model = new Protein(new AtomList("Pisces_REDUCE/"+models[i]), new ResidueExtendedAtoms(ADD_HYDROGENS_AND_FREEZE));
			// Fixing the hydrogen locations
			model.defrost();
			model.freeze(new AtomList.NonHydrogen());
			PutHydrogens.adjustHydrogens(commands, model);
			// Writing the "hydrogenited" PDBs to disk
			try {
				model.atoms().print(new MeshiWriter("Pisces/"+models[i]));
			}
			catch (Exception e) {
				System.out.print("\nThere was a problem writing the output:\n" + e + "\n\nContinueing...\n\n");
			}   */

			// Creating the model
			Protein model = new Protein(new AtomList(models[i]), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
			if (model.atoms().CAFilter().size()>min_Prot_Size) {				
				// Distance matrix and energies
				DistanceMatrix dm = new DistanceMatrix(model.atoms(), cutoffNeighbor+0.5 , 2.0, 4);
				energy = new TotalEnergy(model, dm, energyCreators, commands);
				energy.evaluateAtoms();
				
				// Creating the hydrogen bond list object:
				hbList = new DahiyatHydrogenBondListNoDuplications(dm, model.atoms(), new DahiyatHighAccuracyParamaters());

				// Extracting the solvate data of the protein
				DistanceList dl = dm.nonBondedList();
				for (int c=0 ; c<dl.size(); c++) {
					dis = dl.distanceAt(c).distance();
					if (dis<cutoffNeighbor) {  
						atom1 = dl.distanceAt(c).atom1();
						atom2 = dl.distanceAt(c).atom2();
						if (atom1.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE)==null) 
							atom1.addAttribute(new SolvateExtractionAttribute());
						if (atom2.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE)==null)
							atom2.addAttribute(new SolvateExtractionAttribute());
						attribute1 = (SolvateExtractionAttribute) atom1.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
						attribute2 = (SolvateExtractionAttribute) atom2.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
						hb = hbList.findBondByPolars(atom1, atom2); 
						if (hb!=null) { // Can make SB
							if ((((atom1.type == KNZ) || (atom1.type == RNH) || (atom1.type == TRN)) && 
									((atom2.type == DOD) || (atom2.type == EOE) || (atom2.type == TRO)))      || 
									(((atom2.type == KNZ) || (atom2.type == RNH) || (atom2.type == TRN)) && 
											((atom1.type == DOD) || (atom1.type == EOE) || (atom1.type == TRO)))) { // is SB
								attribute1.addToSumOfSB(hb.hbVal());
								attribute2.addToSumOfSB(hb.hbVal());
							}
							else { // regular HB
								attribute1.addToSumOfHB(hb.hbVal());
								attribute2.addToSumOfHB(hb.hbVal());
							}
						}
						else {  
							if (!atom2.isHydrogen)
								attribute1.addToSumOfNeighbors(1.0);
							if (!atom1.isHydrogen)
								attribute2.addToSumOfNeighbors(1.0);
						}
					}
				}
				
				// Excluding clashing residues
				for (int c=0 ; c<model.atoms().size() ; c++) {
					atom = model.atoms().atomAt(c);
					if (atom.energy()>max_LJ_energy)
						for (int cc=0 ; cc<atom.residue().atoms().size(); cc++)
							atom.residue().atoms().atomAt(cc).addEnergy(100.0);
				}
				

				// Add the protein's data to the entire database
				for (int c=0 ; c<model.atoms().size() ; c++) {
					atom = model.atoms().atomAt(c);
					attribute = (SolvateExtractionAttribute) atom.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
					if (attribute!=null) {
						if (atom.energy()<max_LJ_energy) {
							indNeighbor = findBinIndice(binsNeighbor, attribute.getSumOfNeighbors());
							indHB = findBinIndice(binsHB, attribute.getSumOfHB());
							indSB = findBinIndice(binsSB, attribute.getSumOfSB());
							data[atom.type][indNeighbor][indHB][indSB]++;
						}
						else {
							//System.out.println("Serious clash in: " + atom + "\n");
						}
					}
					else {
						//System.out.println("No Attribute: " + atom + "\n");
					}
				}
			}
		}


		// Outputting
		try{
			DecimalFormat fmt = new DecimalFormat("0.##");    
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
			for (int atomType=0; atomType<Atom.numberOfTypes() ; atomType++) {
				//System.out.println("## type:" + Atom.type(atomType));
				for (int sbc=0; sbc<binsSB.length ; sbc++) {
					//System.out.println("## SBC:" + binsSB[sbc]);
					for (int hbc=0; hbc<binsHB.length ; hbc++) {
						//System.out.print(hbc + "   ");
						for (int cnc=0; cnc<binsNeighbor.length ; cnc++) {
							//System.out.print(data[atomType][cnc][hbc][sbc] + " ");
							bw.write(data[atomType][cnc][hbc][sbc] + " ");
						}
						//System.out.println();
						bw.write("\n");
					}
				}
			}
			bw.close();
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}

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
