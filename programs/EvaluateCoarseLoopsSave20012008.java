package programs;

import java.text.DecimalFormat;
import java.util.StringTokenizer;

import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.ROT1solvation.CBSolvationCreator;
import meshi.energy.ROT1solvation.CentroidSolvationCreator;
import meshi.energy.ROT1solvation.ROT1SolvationCreator;
import meshi.energy.ROT1solvation.ROT1SolvationEnergy;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.compositeTorsions.compositePropensity2DwithPP.CompositePropensity2DCreator;
import meshi.energy.compositeTorsions.ramachandran.RamachandranCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_HighAccuracy_BBonly_Creator;
import meshi.energy.softExcludedVol.SoftExcludedVolCreator;
import meshi.energy.torsionVal.TorsionValCreator;
import meshi.energy.torsionVal.TorsionValEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
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

public class EvaluateCoarseLoopsSave20012008 extends MeshiProgram implements Residues, AtomTypes{ 

	private static CommandList commands; 
	private static String commandsFileName = null;
	private static String modelsFileName = null;  
	private static String dataFileName = null;  
	private static String refFileName = null;  
	private static int resStart = -999;
	private static int resEnd = -999;


	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		DunbrackLib lib = new DunbrackLib(commands, 1.0, 1);
		Protein reference = null;
		Protein model = null;
		DistanceMatrix distanceMatrix = null;
		TotalEnergy energy = null;
		DecimalFormat fmt2 = new DecimalFormat("0.###");


		EnergyCreator[] energyCreatorsFinalEvaluation = {  
				new BondCreator(),
				new AngleCreator(),
				new PlaneCreator(),
				new OutOfPlaneCreator(),
				new SoftExcludedVolCreator(1.0 , 12 , 1.0),
				new SimpleHydrogenBond_Dahiyat_HighAccuracy_BBonly_Creator(1.0),
				new TorsionValCreator(1.0),
				new RamachandranCreator(1.0),
				new CompositePropensity2DCreator(1.0),
				new ROT1SolvationCreator(1.0,"4.0",false),
				new ROT1SolvationCreator(1.0,"4.5",false),
				new ROT1SolvationCreator(1.0,"5.0",false),
				new ROT1SolvationCreator(1.0,"5.5",false),
				new ROT1SolvationCreator(1.0,"6.0",false),
				new ROT1SolvationCreator(1.0,"6.5",false),
				new ROT1SolvationCreator(1.0,"7.0",false),
				new ROT1SolvationCreator(1.0,"8.0",false),
				new ROT1SolvationCreator(1.0,"9.0",false),
				new ROT1SolvationCreator(1.0,"10.0",false),
				new ROT1SolvationCreator(1.0,"12.0",false),
				new ROT1SolvationCreator(1.0,"14.0",false),
				new CBSolvationCreator(1.0,"4.0",false),
				new CBSolvationCreator(1.0,"4.5",false),
				new CBSolvationCreator(1.0,"5.0",false),
				new CBSolvationCreator(1.0,"5.5",false),
				new CBSolvationCreator(1.0,"6.0",false),
				new CBSolvationCreator(1.0,"6.5",false),
				new CBSolvationCreator(1.0,"7.0",false),
				new CBSolvationCreator(1.0,"8.0",false),
				new CBSolvationCreator(1.0,"9.0",false),
				new CBSolvationCreator(1.0,"10.0",false),
				new CBSolvationCreator(1.0,"12.0",false),
				new CBSolvationCreator(1.0,"14.0",false),
				new CentroidSolvationCreator(1.0,"4.0",false),
				new CentroidSolvationCreator(1.0,"4.5",false),
				new CentroidSolvationCreator(1.0,"5.0",false),
				new CentroidSolvationCreator(1.0,"5.5",false),
				new CentroidSolvationCreator(1.0,"6.0",false),
				new CentroidSolvationCreator(1.0,"6.5",false),
				new CentroidSolvationCreator(1.0,"7.0",false),
				new CentroidSolvationCreator(1.0,"8.0",false),
				new CentroidSolvationCreator(1.0,"9.0",false),
				new CentroidSolvationCreator(1.0,"10.0",false),
				new CentroidSolvationCreator(1.0,"12.0",false),
				new CentroidSolvationCreator(1.0,"14.0",false)				
		};	

		 ROT1SolvationEnergy rot1_040 = null;
		 ROT1SolvationEnergy rot1_045 = null;
		 ROT1SolvationEnergy rot1_050 = null;
		 ROT1SolvationEnergy rot1_055 = null;
		 ROT1SolvationEnergy rot1_060 = null;
		 ROT1SolvationEnergy rot1_065 = null;
		 ROT1SolvationEnergy rot1_070 = null;
		 ROT1SolvationEnergy rot1_080 = null;
		 ROT1SolvationEnergy rot1_090 = null;
		 ROT1SolvationEnergy rot1_100 = null;
		 ROT1SolvationEnergy rot1_120 = null;
		 ROT1SolvationEnergy rot1_140 = null;
		 ROT1SolvationEnergy cb_040 = null;
		 ROT1SolvationEnergy cb_045 = null;
		 ROT1SolvationEnergy cb_050 = null;
		 ROT1SolvationEnergy cb_055 = null;
		 ROT1SolvationEnergy cb_060 = null;
		 ROT1SolvationEnergy cb_065 = null;
		 ROT1SolvationEnergy cb_070 = null;
		 ROT1SolvationEnergy cb_080 = null;
		 ROT1SolvationEnergy cb_090 = null;
		 ROT1SolvationEnergy cb_100 = null;
		 ROT1SolvationEnergy cb_120 = null;
		 ROT1SolvationEnergy cb_140 = null;
		 ROT1SolvationEnergy cent_040 = null;
		 ROT1SolvationEnergy cent_045 = null;
		 ROT1SolvationEnergy cent_050 = null;
		 ROT1SolvationEnergy cent_055 = null;
		 ROT1SolvationEnergy cent_060 = null;
		 ROT1SolvationEnergy cent_065 = null;
		 ROT1SolvationEnergy cent_070 = null;
		 ROT1SolvationEnergy cent_080 = null;
		 ROT1SolvationEnergy cent_090 = null;
		 ROT1SolvationEnergy cent_100 = null;
		 ROT1SolvationEnergy cent_120 = null;
		 ROT1SolvationEnergy cent_140 = null;


		// Loading the reference and the background model	
		reference = new Protein((new AtomList(refFileName)).noOXTFilter(), new ResidueExtendedAtoms(ADD_ATOMS));
		RotamericTools.putIntoRot1(reference, new DistanceMatrix(reference.atoms(), 5.5, 2.0, 4), lib);
		model = new Protein((new AtomList(refFileName)).noOXTFilter(), new ResidueExtendedAtoms(ADD_ATOMS));
		for (int cc=0 ; cc<model.atoms().size() ; cc++)
			model.atoms().atomAt(cc).setChain(" ");
		RotamericTools.putIntoRot1(model, new DistanceMatrix(model.atoms(), 5.5, 2.0, 4), lib);
		model.freeze();
		for (int c=resStart ; c<=resEnd ; c++)
			model.residue(c).atoms().defrost();
		distanceMatrix = new DistanceMatrix(model.atoms(), 14.0, 0.1, 4);  
		energy = new TotalEnergy(model, distanceMatrix, energyCreatorsFinalEvaluation, commands);
		TorsionValEnergy torsionTerm1 = (TorsionValEnergy) energy.getEnergyTerm(new TorsionValEnergy());
		rot1_040 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[0]);
		rot1_045 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[1]);
		rot1_050 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[2]);
		rot1_055 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[3]);
		rot1_060 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[4]);
		rot1_065 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[5]);
		rot1_070 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[6]);
		rot1_080 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[7]);
		rot1_090 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[8]);
		rot1_100 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[9]);
		rot1_120 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[10]);
		rot1_140 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[11]);
		cb_040 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[12]);
		cb_045 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[13]);
		cb_050 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[14]);
		cb_055 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[15]);
		cb_060 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[16]);
		cb_065 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[17]);
		cb_070 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[18]);
		cb_080 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[19]);
		cb_090 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[20]);
		cb_100 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[21]);
		cb_120 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[22]);
		cb_140 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[23]);	
		cent_040 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[24]);
		cent_045 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[25]);
		cent_050 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[26]);
		cent_055 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[27]);
		cent_060 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[28]);
		cent_065 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[29]);
		cent_070 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[30]);
		cent_080 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[31]);
		cent_090 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[32]);
		cent_100 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[33]);
		cent_120 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[34]);
		cent_140 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[35]);	

		// Looping on all the models
		String[] models = File2StringArray.f2a(modelsFileName);
		String[] dataFile = File2StringArray.f2a(dataFileName);
		for (int i=-1 ; i<models.length ; i++) {
			try {	
				if (i==-1) {
					// Setting the minimization energy
					String[] phipsiData = File2StringArray.f2a(dataFile[0]);
					for (int c=0 ; c<phipsiData.length ; c++) {
						StringTokenizer st = new StringTokenizer(phipsiData[c]);
						int residueNumber = (new Integer(st.nextToken())).intValue();
						int resType = model.residue(residueNumber).type;
						double phi = (new Double(st.nextToken())).doubleValue();
						double psi = (new Double(st.nextToken())).doubleValue();
						torsionTerm1.getTorsionEnergyElement(residueNumber, "PHI").setTarget(phi);
						torsionTerm1.getTorsionEnergyElement(residueNumber, "PSI").setTarget(psi);
						for (int cc=0 ; cc<lib.getChiMax(resType) ; cc++) {
							double[] rot = lib.getRotamer(resType, phi, psi, 0);
							torsionTerm1.getTorsionEnergyElement(residueNumber, "CHI"+(cc+1)).setTarget(rot[cc]);
						}
					}

					energy.update();
					RotamericTools.putIntoRot1(model, RotamericTools.phipsi(model, distanceMatrix), lib);
					energy.evaluate();
					System.out.print("888888 " + i + " " + calcRMS(model, reference, resStart, resEnd) + " " +
							calcRMSallHeavyAtoms(model, reference, resStart, resEnd) + " " + energy.report(2) + " ");
					
				}
				else {
					// Reading the model
					AtomList loop = new AtomList(models[i]);
					for (int c=0; c<loop.size() ; c++) {
						Atom atom = loop.atomAt(c);
						model.atoms().findAtomInList(atom.name(), atom.residueNumber()).setXYZ(atom.x(), atom.y(), atom.z());
					}

					// Setting the minimization energy
					String[] phipsiData = File2StringArray.f2a(dataFile[i]);
					for (int c=0 ; c<phipsiData.length ; c++) {
						StringTokenizer st = new StringTokenizer(phipsiData[c]);
						int residueNumber = (new Integer(st.nextToken())).intValue();
						int resType = model.residue(residueNumber).type;
						double phi = (new Double(st.nextToken())).doubleValue();
						double psi = (new Double(st.nextToken())).doubleValue();
						torsionTerm1.getTorsionEnergyElement(residueNumber, "PHI").setTarget(phi);
						torsionTerm1.getTorsionEnergyElement(residueNumber, "PSI").setTarget(psi);
						for (int cc=0 ; cc<lib.getChiMax(resType) ; cc++) {
							double[] rot = lib.getRotamer(resType, phi, psi, 0);
							torsionTerm1.getTorsionEnergyElement(residueNumber, "CHI"+(cc+1)).setTarget(rot[cc]);
						}
					}

					energy.update();
					RotamericTools.putIntoRot1(model, RotamericTools.phipsi(model, distanceMatrix), lib);
					energy.evaluate();
					System.out.print("999999 " + i + " " + calcRMS(model, reference, resStart, resEnd) + " " +
							calcRMSallHeavyAtoms(model, reference, resStart, resEnd) + " " + energy.report(2) + " ");
				}
				System.out.print(fmt2.format(rot1_040.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(rot1_045.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(rot1_050.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(rot1_055.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(rot1_060.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(rot1_065.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(rot1_070.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(rot1_080.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(rot1_090.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(rot1_100.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(rot1_120.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(rot1_140.getEnergyHydrophobicWeightedSA()) + "         " + 
						fmt2.format(cb_040.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cb_045.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cb_050.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cb_055.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cb_060.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cb_065.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cb_070.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cb_080.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cb_090.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cb_100.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cb_120.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cb_140.getEnergyHydrophobicWeightedSA()) + "         " + 
						fmt2.format(cent_040.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cent_045.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cent_050.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cent_055.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cent_060.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cent_065.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cent_070.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cent_080.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cent_090.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cent_100.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cent_120.getEnergyHydrophobicWeightedSA()) + " " + 
						fmt2.format(cent_140.getEnergyHydrophobicWeightedSA()) + "         ");
				System.out.print(fmt2.format(rot1_040.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(rot1_045.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(rot1_050.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(rot1_055.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(rot1_060.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(rot1_065.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(rot1_070.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(rot1_080.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(rot1_090.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(rot1_100.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(rot1_120.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(rot1_140.getEnergyPolarWeightedSA()) + "         " + 
						fmt2.format(cb_040.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cb_045.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cb_050.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cb_055.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cb_060.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cb_065.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cb_070.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cb_080.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cb_090.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cb_100.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cb_120.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cb_140.getEnergyPolarWeightedSA()) + "         " + 
						fmt2.format(cent_040.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cent_045.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cent_050.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cent_055.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cent_060.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cent_065.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cent_070.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cent_080.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cent_090.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cent_100.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cent_120.getEnergyPolarWeightedSA()) + " " + 
						fmt2.format(cent_140.getEnergyPolarWeightedSA()) + "         ");
				System.out.println(fmt2.format(rot1_040.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(rot1_045.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(rot1_050.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(rot1_055.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(rot1_060.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(rot1_065.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(rot1_070.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(rot1_080.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(rot1_090.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(rot1_100.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(rot1_120.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(rot1_140.getEnergyRepresentativePolar()) + "         " + 
						fmt2.format(cb_040.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cb_045.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cb_050.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cb_055.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cb_060.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cb_065.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cb_070.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cb_080.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cb_090.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cb_100.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cb_120.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cb_140.getEnergyRepresentativePolar()) + "         " + 
						fmt2.format(cent_040.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cent_045.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cent_050.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cent_055.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cent_060.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cent_065.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cent_070.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cent_080.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cent_090.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cent_100.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cent_120.getEnergyRepresentativePolar()) + " " + 
						fmt2.format(cent_140.getEnergyRepresentativePolar()) + "         ");

			}
			catch (Exception e) {
				System.out.println("Evaluation was not successful.\n");
				e.printStackTrace();
				System.out.println("\nContinueing to next model.\n");
			}
		}

	} // Of main



	private static double calcRMS(Protein prot, Protein ref, int start, int end) {
		double totRms = 0.0;
		for (int c=start; c<=end ; c++) {
			Atom atom = prot.residue(c).atoms().findAtomInList("N",c);
			Atom atomr = ref.residue(c).atoms().findAtomInList("N",c);
			totRms += (atom.x() - atomr.x())*
			(atom.x() - atomr.x()) +
			(atom.y() - atomr.y())*
			(atom.y() - atomr.y()) + 
			(atom.z() - atomr.z())*
			(atom.z() - atomr.z());
			atom = prot.residue(c).atoms().findAtomInList("CA",c);
			atomr = ref.residue(c).atoms().findAtomInList("CA",c);
			totRms += (atom.x() - atomr.x())*
			(atom.x() - atomr.x()) +
			(atom.y() - atomr.y())*
			(atom.y() - atomr.y()) + 
			(atom.z() - atomr.z())*
			(atom.z() - atomr.z());
			atom = prot.residue(c).atoms().findAtomInList("C",c);
			atomr = ref.residue(c).atoms().findAtomInList("C",c);
			totRms += (atom.x() - atomr.x())*
			(atom.x() - atomr.x()) +
			(atom.y() - atomr.y())*
			(atom.y() - atomr.y()) + 
			(atom.z() - atomr.z())*
			(atom.z() - atomr.z());
		}
		return Math.sqrt(totRms/(3*(end-start+1)));
	}

	private static double calcRMSallHeavyAtoms(Protein prot, Protein ref, int start, int end) {
		double totRms = 0.0;
		int ntot = 0;
		for (int c=start; c<=end ; c++) {
			for (int d=0; d<prot.residue(c).atoms().size() ; d++) 
				if (!prot.residue(c).atoms().atomAt(d).isHydrogen) {
					Atom atom = prot.residue(c).atoms().atomAt(d);
					Atom atomr = ref.residue(c).atoms().findAtomInList(prot.residue(c).atoms().atomAt(d).name(),c);
					if (atomr!=null) {
						totRms += (atom.x() - atomr.x())*(atom.x() - atomr.x()) + 
						(atom.y() - atomr.y())*(atom.y() - atomr.y()) +
						(atom.z() - atomr.z())*(atom.z() - atomr.z());
						ntot++;
					}
				}
		}
		return Math.sqrt(totRms/ntot);
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

		dataFileName = getOrderedArgument(args);
		if (dataFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# data file name is "+dataFileName);

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
