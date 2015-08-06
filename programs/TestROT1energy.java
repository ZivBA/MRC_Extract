package programs;

import java.io.IOException;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.rot1Pairwise.CbPairwiseCreator;
import meshi.energy.rot1Pairwise.Rot1PairwiseCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.MinimizerException;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.file.MeshiWriter;
import meshi.util.rotamericTools.RotamericTools;


public class TestROT1energy extends MeshiProgram implements Residues, AtomTypes{ 
    private static CommandList commands; 
   
    private static String commandsFileName = null;
 
    private static String modelFileName = null;  

    private static String loopFileName = null;  

    private static Protein protein;
    
    private static int loopStart=203;

    private static int loopEnd=214;
 
    public static void main(String[] args) throws MinimizerException, LineSearchException{
	init(args); 

	double[] cbEnesNat = new double[loopEnd+1];
	double[] cbEnesLoop = new double[loopEnd+1];
	double[] rot1EnesNat = new double[loopEnd+1];
	double[] rot1EnesLoop = new double[loopEnd+1];
	
	protein = new ExtendedAtomsProtein(modelFileName,DO_NOT_ADD_ATOMS); 
	for (int cc=0 ; cc<protein.atoms().size() ; cc++)
		protein.atoms().atomAt(cc).setChain(" ");
	protein.defrost();

	DunbrackLib lib = new DunbrackLib(commands,1.0,2);
	DistanceMatrix dm = new DistanceMatrix(protein.atoms(), 8.5,  1.0,  4);  
	RotamericTools.putIntoRot1(protein, dm, lib);
	dm = new DistanceMatrix(protein.atoms(), 8.0,  1.0,  4);
	
	try {
		protein.atoms().print(new MeshiWriter(modelFileName+".rot1.pdb"));
	} catch (IOException e1) {
		// TODO Auto-generated catch block
		e1.printStackTrace();
	}

	EnergyCreator[] energyCreators = {  
			new CbPairwiseCreator(1.0)
	};
	TotalEnergy energy = new TotalEnergy(protein, dm, energyCreators, commands);
	energy.evaluate();
	energy.evaluateAtoms();
	for (int res=loopStart; res<=loopEnd ; res++) {
		Atom atom = protein.residue(res).atoms().getAtom("CB");
		if (atom != null)
			cbEnesNat[res] = atom.energy();
	}
	energy.resetAtomEnergies();
	
	AtomList loop = new AtomList(loopFileName);
	for (int c=0; c<loop.size() ; c++) {
		Atom atom = loop.atomAt(c);
		protein.atoms().findAtomInList(atom.name(), atom.residueNumber()).setXYZ(atom.x(), atom.y(), atom.z());
	}
	energy.evaluate();
	RotamericTools.putIntoRot1(protein, dm, lib);
	energy.evaluate();
	energy.evaluateAtoms();
	for (int res=loopStart; res<=loopEnd ; res++) {
		Atom atom = protein.residue(res).atoms().getAtom("CB");
		if (atom != null)
			cbEnesLoop[res] = atom.energy();
	}
	try {
		protein.atoms().print(new MeshiWriter(loopFileName+".rot1.pdb"));
	} catch (IOException e1) {
		// TODO Auto-generated catch block
		e1.printStackTrace();
	}

	
	
	// Now doing the ROT1 rep.
	
	
	
	protein = new ExtendedAtomsProtein(modelFileName,DO_NOT_ADD_ATOMS); 
	for (int cc=0 ; cc<protein.atoms().size() ; cc++)
		protein.atoms().atomAt(cc).setChain(" ");
	protein.defrost();

	dm = new DistanceMatrix(protein.atoms(), 8.5,  1.0,  4);  
	RotamericTools.putIntoRot1(protein, dm, lib);
	for (int res=0 ; res<protein.residues().size() ; res++) {
		if (!protein.residues().residueAt(res).dummy())
			RotamericTools.putCBinSCcenter(protein.residues().residueAt(res));
	}
	dm = new DistanceMatrix(protein.atoms(), 8.0,  1.0,  4);
	
	EnergyCreator[] energyCreators2 = {  
			new Rot1PairwiseCreator(1.0)
	};
	
	energy = new TotalEnergy(protein, dm, energyCreators2, commands);
	energy.evaluate();
	energy.evaluateAtoms();
	for (int res=loopStart; res<=loopEnd ; res++) {
		Atom atom = protein.residue(res).atoms().getAtom("CB");
		if (atom != null)
			rot1EnesNat[res] = atom.energy();
	}
	energy.resetAtomEnergies();
	
	loop = new AtomList(loopFileName);
	for (int c=0; c<loop.size() ; c++) {
		Atom atom = loop.atomAt(c);
		protein.atoms().findAtomInList(atom.name(), atom.residueNumber()).setXYZ(atom.x(), atom.y(), atom.z());
	}
	energy.evaluate();
	RotamericTools.putIntoRot1(protein, dm, lib);
	for (int res=0 ; res<protein.residues().size() ; res++) {
		if (!protein.residues().residueAt(res).dummy())
			RotamericTools.putCBinSCcenter(protein.residues().residueAt(res));
	}
	energy.evaluate();
	System.out.print("*********************************************************");
	energy.evaluateAtoms();
	for (int res=loopStart; res<=loopEnd ; res++) {
		Atom atom = protein.residue(res).atoms().getAtom("CB");
		if (atom != null)
			rot1EnesLoop[res] = atom.energy();
	}
	
	
	System.out.println();
	for (int res=loopStart; res<=loopEnd ; res++) 
		System.out.println(res + " " + (rot1EnesLoop[res]-rot1EnesNat[res]) + " " +
				(cbEnesLoop[res]-cbEnesNat[res]));
	
    }



    /** ================================= init =========================================
     *
     *A static function for parsing of the command line arguments and assigning the 
     *variables commandsFileName, modelFileName and randomNumberSeed with the right inputs. Note that this
     *static method is using parsing functions such as getOrderedArguments that are defined in MeshiProgram
     *that MinimizeProtein inherits.
     **/
     
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
			       "Usage java -Xmx300m SimplestMinimizeProtein <commands file name> <pdb file name> \n"+
			       "                    ******************\n");
			      
	if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));
	commandsFileName = getOrderedArgument(args);
	if (commandsFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# commandsFileName = "+commandsFileName);

	commands = new CommandList(commandsFileName);
	
	modelFileName = getOrderedArgument(args);
	if (modelFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# initial model file name is "+modelFileName);

	loopFileName = getOrderedArgument(args);
	if (loopFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# Loop file name is "+loopFileName);

	initRandom(999);
    }
}
