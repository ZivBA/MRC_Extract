package meshi.solv.parametrization;

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
 * This is a variant on 'ExtractSolvateData', that instead of summarizing on the entire protein,
 * is printing the specific CNC,HBC and SBC values for the atoms in residues: 
 * 110,120,130,140,150,160,170,180 and 190.
 * 
 * @author Nir
 */
public class TestPairwiseEnergies extends MeshiProgram implements Residues, AtomTypes{ 

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

		init(args); 

		// Creators for the energy terms, and separate energy terms
		EnergyCreator[] energyCreators = {
				new LennardJonesCreator(1.0)
		};
		TotalEnergy energy;
		AbstractHydrogenBondList hbList;

		// Going over the models
		String[] models = File2StringArray.f2a(listOfStructures);
		for (int i=0 ; i<models.length ; i++) {
			try {
				System.out.println("Reading: " + prefixName+models[i]+suffixName);

				// Creating the model
				Protein model = new Protein(new AtomList(prefixName+models[i]+suffixName), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
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
								if (!atom.isBackbone) {
									if ((atom.residueNumber()>105) && (atom.residueNumber()<195) && ((atom.residueNumber()%10) == 0)) {
										System.out.println("XXXXX " + models[i] + " " + atom.residueNumber() + " " + atom.residueName() + " " + atom.name + " " +
												attribute.getSumOfNeighbors() + " " + attribute.getSumOfHB() + " " + attribute.getSumOfSB());								
									}
								}
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
				"Usage java -Xmx1000m TestPairwiseEnergies <list of structures> <the distance cutoff of CNC> <Reminder> <prefix of name> <suffix of name>\n"+
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
