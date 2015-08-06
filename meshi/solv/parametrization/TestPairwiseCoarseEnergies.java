package meshi.solv.parametrization;

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
 * This is a variant on 'ExtractCoarseSolvateData', that instead of summarizing on the entire protein,
 * is printing the specific ROT1 - CNC values for the atoms in residues: 
 * 110,120,130,140,150,160,170,180 and 190.
 * 
 * @author Nir
 */

public class TestPairwiseCoarseEnergies extends MeshiProgram implements Residues, AtomTypes{ 

	// We get these user defined parameters from the command line
	// ----------------------------------------------------------
	// The structure database:
	static String listOfStructures = "C:/Users/Nir/Loop_Building_Project/listPISCES.txt";
	// The carbon cutoff (carbon will counted in a sphere with this radius):
	static double cutoffNeighbor = 5.3;    
	// The remider after division by 10 of the residues to take:
	static String prefixName = "";    
	// The remider after division by 10 of the residues to take:
	static String suffixName = "";    


	public static void main(String[] args){
		// ------------------------
		// User defined parameters:
		// ------------------------
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

		// Going over the models
		String[] models = File2StringArray.f2a(listOfStructures);
		for (int i=0 ; i<models.length ; i++) { 
			try {
				System.out.println("Reading: " + prefixName+models[i]+suffixName);

				// Creating the model
				Protein model = null;
				try {
				model = new Protein(new AtomList(prefixName+models[i]+suffixName), new ResidueExtendedAtoms(ADD_ATOMS));
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
//							if (atom2.isBackbone && !atom2.isHydrogen)
							if (atom2.name.equals("CB"))
								attribute1.addToSumOfNeighbors(1.0 - 0.01*dis/cutoffNeighbor);
//							if (atom1.isBackbone && !atom1.isHydrogen)
							if (atom1.name.equals("CB"))
								attribute2.addToSumOfNeighbors(1.0 - 0.01*dis/cutoffNeighbor);
						}
					}

					// Add the protein's data to the entire database
					for (int c=0 ; c<model.atoms().size() ; c++) {
						atom = model.atoms().atomAt(c);
						ty = atom.type;
						attribute = (SolvateExtractionAttribute) atom.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
						if (atom.energy()<max_LJ_energy) {
							if (!atom.isBackbone || atom.name().equals("CB")) {
//								if ((atom.residueNumber()>105) && (atom.residueNumber()<195) && ((atom.residueNumber()%10) == 0)) {
								if (true) {
									System.out.println("XXXXX " + models[i] + " " + atom.residueNumber() + " " + atom.residueName() + " " + atom.name + " " +
											attribute.getSumOfNeighbors());
								}
							}
						}
						else {
							//System.out.println("Serious clash in: " + atom + "\n");
						}
					}


					// ****************************************************
					// ************    ROT1 ROT1 ROT1  ********************
					// ****************************************************
					// Extracting the solvate data for the ROT1
//					RotamericTools.putIntoRot1(model, dm, lib);
//					dm = new DistanceMatrix(model.atoms(), cutoffNeighbor+0.1 , 0.2, 4);
//					for (int c=0; c<model.atoms().size() ; c++)
//						((SolvateExtractionAttribute) model.atoms().atomAt(c).getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE)).reset();		    
//					dl = dm.nonBondedList();
//					for (int c=0 ; c<dl.size(); c++) {
//						dis = dl.distanceAt(c).distance();
//						if ((dis<cutoffNeighbor) && 
//								(Math.abs(dl.distanceAt(c).atom1().residueNumber() - dl.distanceAt(c).atom2().residueNumber())>ignore_local)) {
//							atom1 = dl.distanceAt(c).atom1();
//							atom2 = dl.distanceAt(c).atom2();
//							ty1 = atom1.type;
//							ty2 = atom2.type;
//							attribute1 = (SolvateExtractionAttribute) atom1.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
//							attribute2 = (SolvateExtractionAttribute) atom2.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
//							if (!atom2.isHydrogen)
//								attribute1.addToSumOfNeighbors(1.0 - 0.01*dis/cutoffNeighbor);
//							if (!atom1.isHydrogen)
//								attribute2.addToSumOfNeighbors(1.0 - 0.01*dis/cutoffNeighbor);
//						}
//					}
//
//					// Add the protein's data to the entire database
//					for (int c=0 ; c<model.atoms().size() ; c++) {
//						atom = model.atoms().atomAt(c);
//						ty = atom.type;
//						attribute = (SolvateExtractionAttribute) atom.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
//						if (atom.energy()<max_LJ_energy) {
//							if (!atom.isBackbone || atom.name().equals("CB")) {
////								if ((atom.residueNumber()>105) && (atom.residueNumber()<195) && ((atom.residueNumber()%10) == 0)) {
//								if (true) {
//									System.out.println("YYYYY " + models[i] + " " + atom.residueNumber() + " " + atom.residueName() + " " + atom.name + " " +
//											attribute.getSumOfNeighbors());
//								}
//							}
//						}
//						else {
//							//System.out.println("Serious clash in: " + atom + "\n");
//						}
//					}
				}
			}
			catch (Exception e) {
				System.out.println("Extraction not successful. /n" + e);
			}
		}
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
				"Usage java -Xmx1000m TestPairwiseCoarseEnergies <list of structures> <the distance cutoff of CNC> <Reminder> <prefix of name> <suffix of name>\n"+
		"                    ******************\n");

		listOfStructures = getOrderedArgument(args);
		if (listOfStructures == null) throw new RuntimeException(errorMessage);
		System.out.println("# The structures are taken from: "+listOfStructures);

		String cutoffString = getOrderedArgument(args);
		if (cutoffString== null) throw new RuntimeException(errorMessage);
		cutoffNeighbor = (new Double(cutoffString)).doubleValue();
		System.out.println("# CNC cutoff value: "+ cutoffNeighbor);

		prefixName = getOrderedArgument(args);
		if (prefixName== null) throw new RuntimeException(errorMessage);
		System.out.println("# prefix of name: "+ prefixName);

		suffixName = getOrderedArgument(args);
		if (suffixName== null) throw new RuntimeException(errorMessage);
		System.out.println("# suffix of name: "+ suffixName);

		initRandom(333);
	}
}
