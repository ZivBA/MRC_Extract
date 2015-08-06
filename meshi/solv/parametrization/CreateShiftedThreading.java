package meshi.solv.parametrization;

import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.LennardJones.LennardJones;
import meshi.energy.LennardJones.LennardJonesCreator;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_HighAccuracy_Creator;
import meshi.geometry.DistanceMatrix;
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
import meshi.util.file.MeshiWriter;
import meshi.util.rotamericTools.RotamericTools;
import programs.PutHydrogens;
import programs.SCMOD;

/**
 * This class is a "main" container for the static method:
 * createShiftedThreading(commands, prot, delta)
 * 
 * It will return an AtomList instance of a new protein who's residue N sits on the location of residue
 * N+delta of 'prot'. It is assume that all the atoms are in prot, and that the backbone hydrogens are correct.
 * 
 *
 */

public class CreateShiftedThreading extends MeshiProgram implements Residues, AtomTypes {
	
	private static String commandsFileName = null;
	private static String modelFileName = null;  
	private static String outputFileName = null;  
	private static int delta = 0;  
	
	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args);
		CommandList	commands = new CommandList(commandsFileName);
		Protein prot = new Protein(new AtomList(modelFileName), new ResidueExtendedAtoms(ADD_ATOMS));
		for (int cc=0 ; cc<prot.atoms().size() ; cc++)
			prot.atoms().atomAt(cc).setChain(" ");
		PutHydrogens.adjustHydrogens(commands, prot);
		AtomList shiftedList = createShiftedThreading(commands, prot, delta);
		try {
			shiftedList.print(new MeshiWriter(outputFileName));
		}
		catch (Exception e) {
			System.out.print("\nThere was a problem writing the PDB:\n" + e + "\n\nContinuing...\n\n");
		}		
	}
	
	
	
	public static AtomList createShiftedThreading(CommandList commands, Protein prot, int delta) {
		AtomList newList = new AtomList();
		Protein newProt;
		DistanceMatrix dm;
		DunbrackLib lib  = new DunbrackLib(commands, 1.0 , 90);
		
		// Copying the coordiantes to the new atoms.
		for (int resNum = Math.max(prot.residues().firstNonDummyResidueNumber(), prot.residues().firstNonDummyResidueNumber()-delta) ;
				resNum<= Math.min(((Residue) prot.residues().last()).number,((Residue) prot.residues().last()).number-delta) ;
				resNum++) {
			Atom atom = prot.atoms().findAtomInList("N", resNum+delta);
			if (atom!=null) {
				Atom newAtom = new Atom(atom.x(),atom.y(),atom.z(),"N",
						prot.residue(resNum).nameThreeLetters(),
						resNum,-1);
				newList.add(newAtom);
			}
			atom = prot.atoms().findAtomInList("CA", resNum+delta);
			if (atom!=null) {
				Atom newAtom = new Atom(atom.x(),atom.y(),atom.z(),"CA",
						prot.residue(resNum).nameThreeLetters(),
						resNum,-1);
				newList.add(newAtom);
			}
			atom = prot.atoms().findAtomInList("C", resNum+delta);
			if (atom!=null) {
				Atom newAtom = new Atom(atom.x(),atom.y(),atom.z(),"C",
						prot.residue(resNum).nameThreeLetters(),
						resNum,-1);
				newList.add(newAtom);
			}
			atom = prot.atoms().findAtomInList("O", resNum+delta);
			if (atom!=null) {
				Atom newAtom = new Atom(atom.x(),atom.y(),atom.z(),"O",
						prot.residue(resNum).nameThreeLetters(),
						resNum,-1);
				newList.add(newAtom);
			}
			atom = prot.atoms().findAtomInList("H", resNum+delta);
			if (atom!=null) {
				Atom newAtom = new Atom(atom.x(),atom.y(),atom.z(),"H",
						prot.residue(resNum).nameThreeLetters(),
						resNum,-1);
				newList.add(newAtom);
			}
		}
		
		// Filling the atoms
		newProt = new Protein(newList, new ResidueExtendedAtoms(ADD_ATOMS));

		// adjusting hydrogens and SCMODing
		dm = new DistanceMatrix(newProt.atoms(), 5.5, 2.0, 4);
		double[][] pp = RotamericTools.putIntoRot1(newProt, dm, lib);
		PutHydrogens.adjustHydrogens(commands, prot);
		newProt.defrost();
		try {
			newProt.atoms().print(new MeshiWriter(outputFileName+".rot1.pdb"));
		}
		catch (Exception e) {
			System.out.print("\nThere was a problem writing the PDB:\n" + e + "\n\nContinuing...\n\n");
		}				
		SCMOD.scmod(commands, lib, newProt, 2, 0.0, pp);
		
		// Final minimzation
		EnergyCreator[] energyCreators = {  
				new BondCreator(1.0),
				new AngleCreator(1.0),
				new PlaneCreator(10.0),
				new OutOfPlaneCreator(1.0),
				new LennardJonesCreator(0.4),
				new SimpleHydrogenBond_Dahiyat_HighAccuracy_Creator(3.0),
				new RamachandranSidechainEnergyCreator(1.0)
		};	
		
		newProt.atoms().freeze(new AtomList.BackboneFilter());
		dm = new DistanceMatrix(newProt.atoms(), 5.5, 2.0, 4);  
		TotalEnergy energy = new TotalEnergy(newProt, dm, energyCreators, commands);
		Minimizer minimizer = new LBFGS(energy, 0.001, 50000, 1000);
		try {
			System.out.println(minimizer.minimize());
		}
		catch (Exception e) {
			throw new RuntimeException("ERROR in minimization");
		}		
		
		// printing a LJ values histogram - Start
		// --------------------------------------
		double firstBin = -1.0;
		double lastBin = 4.0;
		double binSpacing = 0.05;
		int[] bins = new int[(int) Math.ceil((lastBin-firstBin)/binSpacing)];
		for (int c=0; c<bins.length ; c++)
			bins[c] = 0;
		energy.resetAtomEnergies();
		LennardJones lj = (LennardJones) energy.getEnergyTerm(new LennardJones());
		lj.evaluateAtoms();
		double Elj = 0.0;
		int ind = 0;
		for (int at=0 ; at<newProt.atoms().size() ; at++) {
			Elj = newProt.atoms().atomAt(at).energy();
			ind = (int) Math.round((Elj-firstBin)/binSpacing);
			if (ind<0)
				ind = 0;
			if (ind>=bins.length)
				ind = bins.length-1;
			bins[ind]++;
		}
		Elj=firstBin;
		for (int c=0; c<bins.length ; c++) {
			System.out.println("000 " + Elj + " " + bins[c]);
			Elj+=binSpacing;
		}
		// --------------------------------------
		// printing the LJ values histogram - END

		return newProt.atoms();
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
				"Usage java -Xmx300m CreateShiftedThreading <commands file name> <model file name> <output file name> <Delta>\n"+
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));
		
		commandsFileName = getOrderedArgument(args);
		if (commandsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# commandsFileName = "+commandsFileName);

		modelFileName = getOrderedArgument(args);
		if (modelFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# model file name is "+modelFileName);

		outputFileName = getOrderedArgument(args);
		if (outputFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# output file name for pdb is "+outputFileName);

		String tmp = getOrderedArgument(args);
		if (tmp== null) throw new RuntimeException(errorMessage);
		delta = (new Integer(tmp)).intValue();
		System.out.println("# delta is:"+delta);
		
		initRandom(999);
	}	
}
