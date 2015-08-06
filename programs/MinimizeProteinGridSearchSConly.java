package programs;

import meshi.applications.prediction.GDTcalculator;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.LennardJones.LennardJonesCreator;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.simpleHPterm.SimpleHP;
import meshi.energy.simpleHPterm.SimpleHPCreator;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBondEnergy;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_Minimization_Creator;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_Minimization_NoDuplications_Creator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.ResidueBuilder;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LBFGS;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.Minimizer;
import meshi.optimizers.MinimizerException;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;
import meshi.util.file.MeshiWriter;
import meshi.util.rotamericTools.RotamericTools;

/**
 *<pre>
 *
 * Unix usage:
 *     java -Xmx300m MinimizeProteinGridSearch <commands file name> <model's list file name> Wbonded Wrg Wlj Whb Whyd Wpol Wramach 
 *
 **/

public class MinimizeProteinGridSearchSConly extends MeshiProgram implements Residues, AtomTypes{ 

	private static CommandList commands; 
	private static String commandsFileName = null;
	private static String modelFileName = null;  
	private static double Wbonded = 0.0;  
	private static double Wrg = 0.0;  
	private static double Wlj = 3.0;  
	private static double Whb = 0.5;  
	private static double Whyd = 1.0;  
	private static double Wpol = 1.0;  
	private static double Wramach = 0.1;  


	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		Protein reference = null;
		Protein model = null;
		DistanceMatrix distanceMatrix = null;
		TotalEnergy energy = null;
		SimpleHydrogenBondEnergy hbTerm,hbTerm1 = null;
		SimpleHP hpTerm = null;
		Minimizer minimizer = null;

		// The creators for the terms
		EnergyCreator[] energyCreators = {  
				new BondCreator(Wbonded),
				new AngleCreator(Wbonded),
				new PlaneCreator(Wrg),
				new OutOfPlaneCreator(Wbonded),
				new LennardJonesCreator(Wlj),
				new SimpleHydrogenBond_Dahiyat_Minimization_NoDuplications_Creator(Whb),
				new SimpleHydrogenBond_Dahiyat_Minimization_Creator(Whb),
				new SimpleHPCreator(Whyd,Wpol,4.25,4.25,true),
				new RamachandranSidechainEnergyCreator(Wramach)
		};	
		
		// The dunbrack lib
		DunbrackLib lib  = new DunbrackLib(commands, 1.0 , 90);

		// Looping on all the models
		String[] proteinNames = File2StringArray.f2a(modelFileName);
		for (int i=0 ; i<proteinNames.length ; i++) {
			// Loading the reference	
			reference = new Protein(new AtomList(proteinNames[i]), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
			double[] rmss = new double[reference.residues().size()];

			// Loading the model	
			model = new Protein(new AtomList(proteinNames[i]), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
			for (int cc=0 ; cc<model.atoms().size() ; cc++)
				model.atoms().atomAt(cc).setChain(" ");			
//			try {
//				model.atoms().print(new MeshiWriter(proteinNames[i]));
//			}
//			catch (Exception e) {
//				System.out.print("\nThere was a problem writing the PDB:\n" + e + "\n\nContinuing...\n\n");
//			}
			distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
			double[][] pp = RotamericTools.phipsi(model, distanceMatrix);
			for (int c=0 ; c<pp.length ; c++) 
				if (pp[c]!=null)
					if ((model.residue(c).type!=ALA) && (model.residue(c).type!=GLY)) {
						double[] rotamer = RotamericTools.nearestRot(lib, model.residue(c), pp[c][0], pp[c][1]);
						if (rotamer!=null) {
							ResidueBuilder.build(model.residue(c), model.residue(c).type, rotamer);
						}
						else {
							model.residue(c).freeze();
						}
					}
			model.freeze(new AtomList.BackboneFilter());
			for (int rr=0 ; rr<reference.residues().size() ; rr++) {
				if (!reference.residues().residueAt(rr).dummy() && (reference.residues().residueAt(rr).type!=ALA) &&
						(reference.residues().residueAt(rr).type!=GLY)) {
					try {
						rmss[rr] = RotamericTools.calcRMS(reference.residues().residueAt(rr), model.residue(reference.residues().residueAt(rr).number));
					}
					catch (Exception e) {
						// Do nothing
					}
				}
			}
			try {
				model.atoms().print(new MeshiWriter(proteinNames[i]+".rot.pdb"));
			}
			catch (Exception e) {
				System.out.print("\nThere was a problem writing the PDB:\n" + e + "\n\nContinuing...\n\n");
			}
			
			String header = "999999" + i + "9999 ";
			System.out.println(header + i + " 0000 " + i);
			System.out.println(header + i + " 1111 " + reference.atoms().CAFilter().getRms(model.atoms().CAFilter()) + 
					" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),0.5,1.0,2.0,4.0) +	
					" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),1.0,2.0,4.0,8.0) + " " +	
					reference.atoms().noOXTFilter().filter(new AtomList.NonHydrogen()).getRms(model.atoms().noOXTFilter().filter(new AtomList.NonHydrogen())));
			/* Initial energies of the model */ 
			distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
			energy = new TotalEnergy(model, distanceMatrix, energyCreators, commands);
			energy.evaluate();
			hbTerm = (SimpleHydrogenBondEnergy) energy.getEnergyTerms(new SimpleHydrogenBondEnergy())[0];
			hbTerm1 = (SimpleHydrogenBondEnergy) energy.getEnergyTerms(new SimpleHydrogenBondEnergy())[1];
			hbTerm.on();
			hbTerm1.off();
			hpTerm = (SimpleHP) energy.getEnergyTerm(new SimpleHP());
			System.out.println(header + i + " 2222 " + energy.report(2) + " " + 
					hpTerm.hydrophobicEnergy() + " " + hpTerm.hydrophilicEnergy() + " " + hbTerm.hbEnergy());

			// Energy and RMS - After minimization
			minimizer = new LBFGS(energy, 0.001, 275, 1000);
			try {
				System.out.println(minimizer.minimize());
			}
			catch (Exception e) {
				System.out.println("ERROR in minimization");
			}
			hbTerm.off();
			hbTerm1.on();
			minimizer = new LBFGS(energy, 0.001, 50000, 1000);
			try {
				System.out.println(minimizer.minimize());
			}
			catch (Exception e) {
				System.out.println("ERROR in minimization");
			}			
			System.out.println(header + i + " 3333 " + reference.atoms().CAFilter().getRms(model.atoms().CAFilter()) + 
					" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),0.5,1.0,2.0,4.0) +	
					" " + GDTcalculator.gdt(reference.atoms(),model.atoms(),1.0,2.0,4.0,8.0) + " " +	
					reference.atoms().noOXTFilter().filter(new AtomList.NonHydrogen()).getRms(model.atoms().noOXTFilter().filter(new AtomList.NonHydrogen())));
			/* Initial energies of the model */ 
			hbTerm = (SimpleHydrogenBondEnergy) energy.getEnergyTerms(new SimpleHydrogenBondEnergy())[0];
			hbTerm1 = (SimpleHydrogenBondEnergy) energy.getEnergyTerms(new SimpleHydrogenBondEnergy())[1];
			hpTerm = (SimpleHP) energy.getEnergyTerm(new SimpleHP());
			System.out.println(header + i + " 4444 " + energy.report(2) + " " + 
					hpTerm.hydrophobicEnergy() + " " + hpTerm.hydrophilicEnergy() + " " + hbTerm1.hbEnergy());
			for (int rr=0 ; rr<reference.residues().size() ; rr++) {
				if (!reference.residues().residueAt(rr).dummy() && (reference.residues().residueAt(rr).type!=ALA) &&
						(reference.residues().residueAt(rr).type!=GLY)) {
					try {
						rmss[rr] -= RotamericTools.calcRMS(reference.residues().residueAt(rr), model.residue(reference.residues().residueAt(rr).number));
						System.out.println(reference.residues().residueAt(rr) + ":  " + rmss[rr]);
					}
					catch (Exception e) {
						// Do nothing
					}
				}
			}
			try {
				model.atoms().print(new MeshiWriter(proteinNames[i]+".min.pdb"));
			}
			catch (Exception e) {
				System.out.print("\nThere was a problem writing the PDB:\n" + e + "\n\nContinuing...\n\n");
			}
		}

	} // Of main



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
				"Usage java -Xmx300m MinimizeProteinGridSearch <commands file name> <model's list file name> Wbonded Wrg Wlj Whb Whyd Wpol Wramach\n"+
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));
		commandsFileName = getOrderedArgument(args);
		if (commandsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# commandsFileName = "+commandsFileName);

		commands = new CommandList(commandsFileName);

		modelFileName = getOrderedArgument(args);
		if (modelFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# loop model file name is "+modelFileName);

		String tmp = getOrderedArgument(args);
		if (tmp== null) throw new RuntimeException(errorMessage);
		Wbonded = (new Double(tmp)).doubleValue();
		System.out.println("# Wbonded is:"+Wbonded);

		tmp = getOrderedArgument(args);
		if (tmp== null) throw new RuntimeException(errorMessage);
		Wrg = (new Double(tmp)).doubleValue();
		System.out.println("# Wrg is:"+Wrg);

		tmp = getOrderedArgument(args);
		if (tmp== null) throw new RuntimeException(errorMessage);
		Wlj = (new Double(tmp)).doubleValue();
		System.out.println("# Wlj is:"+Wlj);

		tmp = getOrderedArgument(args);
		if (tmp== null) throw new RuntimeException(errorMessage);
		Whb = (new Double(tmp)).doubleValue();
		System.out.println("# Whb is:"+Whb);

		tmp = getOrderedArgument(args);
		if (tmp== null) throw new RuntimeException(errorMessage);
		Whyd = (new Double(tmp)).doubleValue();
		System.out.println("# Whyd is:"+Whyd);

		tmp = getOrderedArgument(args);
		if (tmp== null) throw new RuntimeException(errorMessage);
		Wpol = (new Double(tmp)).doubleValue();
		System.out.println("# Wpol is:"+Wpol);

		tmp = getOrderedArgument(args);
		if (tmp== null) throw new RuntimeException(errorMessage);
		Wramach = (new Double(tmp)).doubleValue();
		System.out.println("# Wramach is:"+Wramach);

		initRandom(999);
	}
}
