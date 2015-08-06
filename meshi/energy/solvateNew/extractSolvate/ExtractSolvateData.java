package meshi.energy.solvateNew.extractSolvate;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.DecimalFormat;

import meshi.energy.excludedVolumeImprovedDistance.EVenergyParametersList;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.hydrogenBonds.AbstractHydrogenBond;
import meshi.molecularElements.hydrogenBonds.AbstractHydrogenBondList;
import meshi.molecularElements.hydrogenBonds.newParameters.HydrogenBondDahiyatList;
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
 * 1. Number of carbon neighbors.
 * 2. Number of hydrogen bonds.
 * 3. Number of salt bridges.
 * 
 * @author Nir
 */
public class ExtractSolvateData extends MeshiProgram implements Residues, AtomTypes{ 
 
    public static void main(String[] args){
    	
    	// ------------------------
    	// User defined parameters:
    	// ------------------------
    	// The structure database:
    	String listOfStructures = "C:/Users/Nir/Loop_Building_Project/listPISCES.txt";
    	// The carbon bins:
    	double[] binsNeighbor = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100};
    	// The hb bins:
    	double[] binsHB = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    	//double[] binsHB = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};
    	// The salt bridge bins:
    	double[] binsSB = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    	//double[] binsSB = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};
    	// The carbon cutoff (carbon will counted in a sphere with this radius):
    	double cutoffNeighbor = 5.3;    	
    	// Minimum size protein:
    	int min_Prot_Size = 69;    	
    	// Output file:
    	String outputFile = "All_Atom_Solvate_Data_5.3"; 
    	// EV parameter file:
    	String evParamasFile = "C:/Users/Nir/Nir_Programs_eclipse/MESHI_E/meshicurrent/parameters/meshiPotential/improvedEVparameters.dat";
    	// CommandList:
    	CommandList commands = new CommandList("commands");
    	// ----------------------------------------
    	
    	// Local variables
		double dis;
		int ty,ty1,ty2;
		Atom atom,atom1,atom2;
		SolvateExtractionAttribute attribute,attribute1,attribute2;
		AbstractHydrogenBond hb;
		int indNeighbor,indHB,indSB;
    	
    	init(args); 
    	
    	// Creating the EV parameters for the Hb list object
    	EVenergyParametersList evParams = new EVenergyParametersList(evParamasFile);
    	
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
			Protein model = new Protein(new AtomList("Pisces/"+models[i]), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
			DistanceMatrix dm = new DistanceMatrix(model.atoms(), cutoffNeighbor+0.5 , 2.0, 4);
			if (model.atoms().CAFilter().size()>min_Prot_Size) {				
		    	// Creating the hydrogen bond list object:
		    	AbstractHydrogenBondList hbList = new HydrogenBondDahiyatList(dm,model.atoms(),evParams);

		    	// Extracting the solvate data of the protein
				DistanceList dl = dm.nonBondedList();
				for (int c=0 ; c<dl.size(); c++) {
					dis = dl.distanceAt(c).distance();
					if ((dis<cutoffNeighbor) && 
							(dl.distanceAt(c).atom1().residueNumber() != dl.distanceAt(c).atom2().residueNumber())) {
						atom1 = dl.distanceAt(c).atom1();
						atom2 = dl.distanceAt(c).atom2();
						ty1 = atom1.type;
						ty2 = atom2.type;
						if (atom1.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE)==null) 
							atom1.addAttribute(new SolvateExtractionAttribute());
						if (atom2.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE)==null)
							atom2.addAttribute(new SolvateExtractionAttribute());
						attribute1 = (SolvateExtractionAttribute) atom1.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
						attribute2 = (SolvateExtractionAttribute) atom2.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
						hb = hbList.findBondByPolars(atom1, atom2); 
						if (hb!=null) { // Can make HB or SB
							if ((((atom1.type == KNZ) || (atom1.type == RNH) || (atom1.type == TRN)) && 
						    	     ((atom2.type == DOD) || (atom2.type == EOE) || (atom2.type == TRO)))      || 
						    	     (((atom2.type == KNZ) || (atom2.type == RNH) || (atom2.type == TRN)) && 
						    	     ((atom1.type == DOD) || (atom1.type == EOE) || (atom1.type == TRO)))) { // is SB
								attribute1.addToSumOfSB(hb.hbVal());
								attribute2.addToSumOfSB(hb.hbVal());
							}
							else { // is HB
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
				
				// Add the protein's data to the entire database
				for (int c=0 ; c<model.atoms().size() ; c++) {
					atom = model.atoms().atomAt(c);
					ty = atom.type;
					attribute = (SolvateExtractionAttribute) atom.getAttribute(MeshiAttribute.SOLVATE_EXTRACTION_ATTRIBUTE);
					if (attribute==null) {
						//System.out.println("No Attribute: " + atom + "\n");
					}
					else {
						indNeighbor = findBinIndice(binsNeighbor, attribute.getSumOfNeighbors());
						indHB = findBinIndice(binsHB, attribute.getSumOfHB());
						indSB = findBinIndice(binsSB, attribute.getSumOfSB());
						data[ty][indNeighbor][indHB][indSB]++;
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

	initRandom(333);
    }
}
