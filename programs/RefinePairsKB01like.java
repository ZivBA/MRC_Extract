package programs;

import java.io.File;
import java.util.StringTokenizer;

import meshi.applications.corpus.Corpus;
//import meshi.applications.corpus.CorrectRottenRamachandrans;
import meshi.applications.prediction.GDTcalculator;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.LennardJones.LennardJonesCreator;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.compositeTorsions.ramachandran.RamachandranCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBondEnergy;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_Minimization_Creator;
import meshi.energy.tether.TetherCreator;
import meshi.energy.torsionVal.TorsionValCreator;
import meshi.energy.torsionVal.TorsionValEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.ResidueBuilder;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LBFGS;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.Minimizer;
import meshi.optimizers.MinimizerException;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.external.ComplexMESHIconversion;
import meshi.util.file.File2StringArray;
import meshi.util.file.MeshiWriter;
import meshi.util.rotamericTools.RotamericTools;

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

public class RefinePairsKB01like extends MeshiProgram implements Residues, AtomTypes{ 

	private static CommandList commands; 
	private static String commandsFileName = null;
	private static String corpusFileName = null;	
	private static String modelsFileName = null;  
	private static String refFileName = null;  
	private static String outputFileName = null;  
	private static String letters = null;


	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		Protein reference = null;
		Protein model = null;
		DistanceMatrix distanceMatrix = null;
		TotalEnergy energy = null;
		SimpleHydrogenBondEnergy hbTerm = null;
		Minimizer minimizer = null;		
		boolean fullMatchRefModel = false;

//		Corpus corpus = new Corpus(corpusFileName);
		DunbrackLib lib  = new DunbrackLib(commands, 1.0 , 100);

		// The creators for the terms
		EnergyCreator[] energyCreators1 = {  
				new BondCreator(1.0),
				new AngleCreator(1.0),
				new PlaneCreator(10.0),
				new OutOfPlaneCreator(1.0),
				new LennardJonesCreator(0.4),
				new SimpleHydrogenBond_Dahiyat_Minimization_Creator(1.0),
				//new RamachandranSidechainEnergyCreator(1.0),
				new TorsionValCreator(20),
				new RamachandranCreator(1.0),
				new TetherCreator(0.5, new AtomList.ClassCAFilter())  //.ClassCaCbFilter()) 
		};	


		// Looping on all the models
		String[] models = File2StringArray.f2a(modelsFileName);
		String[] refs = File2StringArray.f2a(refFileName);
		String[] outputs = File2StringArray.f2a(outputFileName);
		String[] lettersA = File2StringArray.f2a(letters);
		for (int i=0 ; i<models.length ; i++) { 		
			// Reading the model
			System.out.println("Trying to minimize: " + models[i]);
			Atom.resetNumberOfAtoms();
			try {
				String chainN1 = ""+lettersA[i].charAt(0);
				String chainN2 = ""+lettersA[i].charAt(2);
				// Loading the reference	
				AtomList refList = new AtomList(refs[i]);
				refList.chainFilter(chainN1).setChain("X");
				refList.chainFilter(chainN2).setChain("Y");
				refList.chainFilter("X").setChain("A");
				refList.chainFilter("Y").setChain("B");
				reference = ComplexMESHIconversion.complex2meshi(refList);
				// Loading the initial model
				AtomList modelList = new AtomList(models[i]);
				modelList.chainFilter(chainN1).setChain("X");
				modelList.chainFilter(chainN2).setChain("Y");
				modelList.chainFilter("X").setChain("A");
				modelList.chainFilter("Y").setChain("B");
				model = ComplexMESHIconversion.complex2meshi(modelList);
				if (getMatchingAtoms(reference, model.atoms())==null)
					fullMatchRefModel = false;
				else
					fullMatchRefModel = true;
				for (int cc=0 ; cc<model.atoms().size() ; cc++)
					model.atoms().atomAt(cc).setChain(" ");
				PutHydrogens.adjustHydrogens(commands, model);
				model.defrost();


				// Energy and RMS - Before minimization
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
				distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
				energy = new TotalEnergy(model, distanceMatrix, energyCreators1, commands);
				energy.evaluate();
				hbTerm = (SimpleHydrogenBondEnergy) energy.getEnergyTerm(new SimpleHydrogenBondEnergy());
				System.out.println("999999 " + i + " 2222 " + energy.report(2) + " " + 
						hbTerm.hbEnergy());					

				// Correcting the bad phi-psi
				if ((new File(models[i]+".noRottenMinimialSCMOD")).exists()) {
					System.out.print("\nSCMOD already run. Continuing...\n\n");
					model = new Protein(new AtomList(models[i]+".noRottenMinimialSCMOD"), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
				}
				else {
					double[][] pp = RotamericTools.phipsi(model, distanceMatrix);
//					CorrectRottenRamachandrans corr = new CorrectRottenRamachandrans(commands , corpus, lib);
//					double[][] noRot = corr.detectAndCorrectRottenResidues(model);
//					// Fixing the changed side-chains
//					model.freeze();
//					for (int c=0 ; c<noRot.length; c++) {
//						pp[(int) noRot[c][0]][0] = noRot[c][1]; // new Phi
//						pp[(int) noRot[c][0]][1] = noRot[c][2]; // new Psi
//						model.residue((int) noRot[c][0]).atoms().defrost();
//					}
//					//model.defrost();  Uncomment this for FullSCMOD
//					SCMOD.scmod(commands, lib, model, 2, 0.0, pp);
					model.defrost();
					for (int res=0 ; res<model.residues().size() ; res++) {
						Residue residue = model.residues().residueAt(res);
						if ((residue.type>0) && (residue.type<20) && (residue.type!=GLY) && (residue.type!=ALA)) {
							double[] nearestRot = RotamericTools.nearestRot(lib, residue, pp[residue.number][0] , pp[residue.number][1]);
							if (nearestRot!=null) {
								ResidueBuilder.build(residue, residue.type,  nearestRot);
							}
						}
					}

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
					distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
					energy = new TotalEnergy(model, distanceMatrix, energyCreators1, commands);
					energy.evaluate();
					hbTerm = (SimpleHydrogenBondEnergy) energy.getEnergyTerm(new SimpleHydrogenBondEnergy());
					System.out.println("999999 " + i + " 4444 " + energy.report(2) + " " + 
							hbTerm.hbEnergy());					
				}

				// Straightening up after rotten correction
				distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);
				energy = new TotalEnergy(model, distanceMatrix, energyCreators1, commands);
				((TorsionValEnergy) (energy.getEnergyTerm(new TorsionValEnergy()))).setChisTargetsAutomaticly();
				energy.evaluate();
				minimizer = new LBFGS(energy, 0.05, 20000, 100);
				System.out.println(minimizer.minimize());
				if (fullMatchRefModel)
					System.out.println("999999 " + i + " 5555 " + reference.atoms().CAFilter().getRms(getMatchingAtoms(reference,model.atoms()).CAFilter()) + 
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),0.5,1.0,2.0,4.0) +	
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),1.0,2.0,4.0,8.0) + " " +	
							reference.atoms().noOXTFilter().filter(new AtomList.NonHydrogen()).getRms(getMatchingAtoms(reference,model.atoms()).noOXTFilter().filter(new AtomList.NonHydrogen()))
							+ " " + getPercentageOfMatchingAtoms(reference.atoms(), model.atoms()));
				else
					System.out.println("999999 " + i + " 5555 " + (-1.0) + 
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),0.5,1.0,2.0,4.0) +	
							" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),1.0,2.0,4.0,8.0) + " " +	
							(-1.0)
							+ " " + getPercentageOfMatchingAtoms(reference.atoms(), model.atoms()));
				distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
				energy = new TotalEnergy(model, distanceMatrix, energyCreators1, commands);
				energy.evaluate();
				hbTerm = (SimpleHydrogenBondEnergy) energy.getEnergyTerm(new SimpleHydrogenBondEnergy());
				System.out.println("999999 " + i + " 6666 " + energy.report(2) + " " + 
						hbTerm.hbEnergy());											
				try {
					AtomList outModel = ComplexMESHIconversion.MEHSI2complex(model);
					outModel.chainFilter("A").setChain("X");
					outModel.chainFilter("B").setChain("Y");
					outModel.chainFilter("X").setChain(chainN1);
					outModel.chainFilter("Y").setChain(chainN2);
					outModel.noOXTFilter().filter(new AtomList.NonHydrogen()).print(new MeshiWriter(outputs[i]));
				}
				catch (Exception e) {
					System.out.print("\nThere was a problem writing the minimized output:\n" + e + "\n\nContinuing...\n\n");
				}
			} 
			catch (Exception e) {
				System.out.println("SKIPPING: minimization was not successful: " + models[i]);
			}
		} // Loop on all the models
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


		String errorMessage = ("\n                  ******************\n"+
				"Usage java -Xmx300m RefineForKB01 <commands file name> <frag corpus file name> <file with list of initial pdbs> <file with list of reference pdbs> <file with list of output names>\n" + 
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));
		commandsFileName = getOrderedArgument(args);
		if (commandsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# commandsFileName = "+commandsFileName);

		commands = new CommandList(commandsFileName);

		corpusFileName = getOrderedArgument(args);
		if (corpusFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Frag corpus file name is "+corpusFileName);

		modelsFileName = getOrderedArgument(args);
		if (modelsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# initial model file name is "+modelsFileName);

		refFileName = getOrderedArgument(args);
		if (refFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# reference file name is "+refFileName);

		outputFileName = getOrderedArgument(args);
		if (outputFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# output file name is "+outputFileName);

		letters = getOrderedArgument(args);
		if (letters == null) throw new RuntimeException(errorMessage);
		System.out.println("# letters file name is "+letters);
		
		initRandom(999);
	}
}
