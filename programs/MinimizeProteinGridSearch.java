package programs;

import meshi.applications.prediction.GDTcalculator;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.LennardJones.LennardJonesCreator;
import meshi.energy.angle.AngleCreator;
import meshi.energy.angle.AngleEnergy;
import meshi.energy.bond.BondCreator;
import meshi.energy.bond.BondEnergy;
import meshi.energy.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.linearRG.LinearRgCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.outOfPlane.OutOfPlaneEnergy;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.plane.PlaneEnergy;
import meshi.energy.simpleHPterm.SimpleHP;
import meshi.energy.simpleHPterm.SimpleHPCreator;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBondEnergy;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_HighAccuracy_Creator;
import meshi.energy.tether.TetherCreator;
import meshi.energy.tether.TetherEnergy;
import meshi.energy.torsionSpaceMinimization.TotalEnergyTorsionSpace;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LBFGS;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.Minimizer;
import meshi.optimizers.MinimizerException;
import meshi.optimizers.SteepestDecent;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;
import meshi.util.filters.KolDichfin;

/**
 *<pre>
 *
 * Unix usage:
 *     java -Xmx300m MinimizeProteinGridSearch <commands file name> <model's list file name> Wbonded Wrg Wlj Whb Whyd Wpol Wramach 
 *
 **/

public class MinimizeProteinGridSearch extends MeshiProgram implements Residues, AtomTypes{ 

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
		SimpleHydrogenBondEnergy hbTerm = null;
		SimpleHP hpTerm = null;
		TetherEnergy tetherTerm = null;
		Minimizer minimizer = null;

		// The creators for the terms
		EnergyCreator[] energyCreators = {  
				new BondCreator(Wbonded),
				new AngleCreator(Wbonded),
				new PlaneCreator(Wbonded),
				new OutOfPlaneCreator(Wbonded),
				new LinearRgCreator(Wrg),
				new LennardJonesCreator(Wlj),
				new SimpleHydrogenBond_Dahiyat_HighAccuracy_Creator(Whb),
				new SimpleHPCreator(Whyd,Wpol,4.25,4.25,true),
				new RamachandranSidechainEnergyCreator(Wramach),
				new TetherCreator(5.0,new KolDichfin())
		};	

		// Looping on all the models
		String[] proteinNames = File2StringArray.f2a(modelFileName);
		for (int i=0 ; i<proteinNames.length ; i++) {
			// Loading the reference	
			reference = new Protein(new AtomList(proteinNames[i]), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));

			// Loading the model	
			model = new Protein(new AtomList(proteinNames[i]), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
			
			String header = "999999" + i + "9999 ";
			System.out.println(header + i + " 0000 " + i);
			System.out.println(header + i + " 1111 0.0 1.0 1.0 0.0");
			/* Initial energies of the model */ 
			distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
			energy = new TotalEnergyTorsionSpace(model, distanceMatrix, energyCreators, commands);
			energy.evaluate();
			hbTerm = (SimpleHydrogenBondEnergy) energy.getEnergyTerm(new SimpleHydrogenBondEnergy());
			hpTerm = (SimpleHP) energy.getEnergyTerm(new SimpleHP());
			tetherTerm = (TetherEnergy) energy.getEnergyTerm(new TetherEnergy());
			System.out.println(header + i + " 2222 " + energy.report(2) + " " + 
					hpTerm.hydrophobicEnergy() + " " + hpTerm.hydrophilicEnergy() + " " + hbTerm.hbEnergy());

			energy.getEnergyTerm(new BondEnergy()).off();
			energy.getEnergyTerm(new AngleEnergy()).off();
			energy.getEnergyTerm(new PlaneEnergy()).off();
			energy.getEnergyTerm(new OutOfPlaneEnergy()).off();
//			energy.resetAtomEnergies();
//			energy.evaluateAtoms();

			// Energy and RMS - After minimization
			minimizer = new SteepestDecent(energy, 0.001, 200, 100);
			try {
				System.out.println(minimizer.minimize());
			}
			catch (Exception e) {
				System.out.println("ERROR in minimization");
			}
			tetherTerm.off();
			minimizer = new LBFGS(energy, 0.001, 10000, 1000);
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
			hbTerm = (SimpleHydrogenBondEnergy) energy.getEnergyTerm(new SimpleHydrogenBondEnergy());
			hpTerm = (SimpleHP) energy.getEnergyTerm(new SimpleHP());

			energy.getEnergyTerm(new BondEnergy()).on();
			energy.getEnergyTerm(new AngleEnergy()).on();
			energy.getEnergyTerm(new PlaneEnergy()).on();
			energy.getEnergyTerm(new OutOfPlaneEnergy()).on();
			energy.evaluate();
			
			System.out.println(header + i + " 4444 " + energy.report(2) + " " + 
					hpTerm.hydrophobicEnergy() + " " + hpTerm.hydrophilicEnergy() + " " + hbTerm.hbEnergy());
//			try {
//				model.atoms().print(new MeshiWriter(proteinNames[i]+".min.pdb"));
//			}
//			catch (Exception e) {
//				System.out.print("\nThere was a problem writing the PDB:\n" + e + "\n\nContinuing...\n\n");
//			}
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
