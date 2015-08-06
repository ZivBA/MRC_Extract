package programs;

import java.io.IOException;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.applications.prediction.GDTcalculator;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.compositeTorsions.compositePropensity2DwithPP.CompositePropensity2DCreator;
import meshi.energy.compositeTorsions.compositePropensity2DwithPP.CompositePropensity2DEnergy;
import meshi.energy.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.inflate.InflateCreator;
import meshi.energy.linearRG.LinearRgCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.softExcludedVol.SoftExcludedVolCreator;
import meshi.energy.solvate.SolvateCreatorHBforMinimization;
import meshi.energy.solvate.SolvateEnergy;
import meshi.energy.tether.TetherCreator;
import meshi.energy.tether.TetherEnergy;
import meshi.geometry.DistanceMatrix;
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


public class MinimizeProteinInflate extends MeshiProgram implements Residues, 
							     AtomTypes{ /**
									 * The implemented
									 * interfaces defines the 
									 * names of atom and residue 
									 * types. 
									 **/
    /**
     * Reads, parse and stores the contents of the commands file. 
     **/ 
    private static CommandList commands; 
   
    /**
     * A string with the name of the command file.
     **/
    private static String commandsFileName = null;
 
    /**
     * A string with the name of the pdb file to minimize.
     **/
    private static String modelFileName = null;  

    /**
     * A string with the name of the pdb file of reference.
     **/
    private static String refFileName = null;  


    private static int modelNum = 1;  
    private static int numOfSCMODiterations = 0;  
    private static int withSCMOD = 0;  
    private static double randomFactorSCMOD = 0.0;  
    private static double Wrg = 0.0;  
    private static double Wev = 3.0;  
    private static double Wsolv = 0.5;  
    private static double Whb = 1.0;  
    private static double Wprop = 1.0;  
    private static double Wramach = 0.1;  
    private static double Wsc = 0.0;  
	private static double Winf = 0.0;
	private static double infTar = 0.0;


    /**
     * The minimized protein object.
     **/
    private static Protein origProt;
    
    private static Protein protein;

    private static Protein reference;  
 
    private static Protein copyReference;  

    public static void main(String[] args) throws IOException,MinimizerException, LineSearchException{
	/**
	 * A static function for parsing of the command line arguments and assigning the 
	 * variables commandsFileName, modelFileName and randomNumberSeed with the right values.
	 * Note that this method is using parsing functions such as getOrderedArguments that are 
	 * defined in MeshiProgram that MinimizeProtein extends.
	 **/
	init(args); 

	AtomList chosenList = null;
	int iteration = 0;
	int bestList = -1;
	double bestEnergy = 1e10;
	double tmpEnergy = 0;
	double LocWprop = 0.5;
	double LocWsolv = 0.3;
	double LocWhb = 0.0;
	reference = new Protein(refFileName);
	copyReference = new ExtendedAtomsProtein(refFileName,ADD_ATOMS_AND_FREEZE); 	
	origProt = new Protein(getMatchingAtoms(reference.atoms().noOXTFilter(),
		new AtomList(modelFileName)) , new ResidueExtendedAtoms(ADD_ATOMS_AND_FREEZE));
	for (int cc=0 ; cc<origProt.atoms().size() ; cc++)
		origProt.atoms().atomAt(cc).setChain("A");
	


	EnergyCreator[] energyCreators = {  
	    new BondCreator(),
	    new AngleCreator(),
	    new PlaneCreator(),
	    new OutOfPlaneCreator(),
	    //new ExcludedVolCreator(),
		//new SolvateNoHBCreator(0.0),
		new LinearRgCreator(Wrg),
		new SoftExcludedVolCreator(Wev , 3 , 1.0),
		//new SolvateCreatorWideAttraction(Wsolv,Whb,95.0,100.0,85.0,90.0),
		new SolvateCreatorHBforMinimization(Wsolv,Wsolv,Wsolv,Whb),
		new CompositePropensity2DCreator(Wprop),
		new RamachandranSidechainEnergyCreator(Wramach),
		new InflateCreator(Winf, infTar)			    
	};

	EnergyCreator[] energyCreatorsNoInf = {  
	    new BondCreator(),
	    new AngleCreator(),
	    new PlaneCreator(),
	    new OutOfPlaneCreator(),
	    //new ExcludedVolCreator(),
	    new TetherCreator(1.0 , new AtomList.ClassCAFilter()),
		//new SolvateNoHBCreator(0.0),
		new LinearRgCreator(Wrg),
		new SoftExcludedVolCreator(Wev , 3 , 1.0),
		//new SolvateCreatorWideAttraction(Wsolv,Whb,95.0,100.0,85.0,90.0),
		new SolvateCreatorHBforMinimization(Wsolv,Wsolv,Wsolv,Whb),
		new CompositePropensity2DCreator(Wprop),
		new RamachandranSidechainEnergyCreator(Wramach),
	};	

	System.out.println("************* Phase I - part of the protein is frozen *********");
	System.out.println("************* Minimizing the hydrogens of model *********");
	DistanceMatrix distanceMatrix = new DistanceMatrix(origProt.atoms(), 5.5, 2.0,4); 
	TotalEnergy energy = new TotalEnergy(origProt, distanceMatrix, energyCreatorsNoInf, commands);
	Minimizer minimizer = new LBFGS(energy, commands);  
	minimizer.minimize();
	origProt.defrost();
	if (withSCMOD>0) {
	    DunbrackLib lib = new DunbrackLib(commands,1.0,30); 
    	SCMOD.scmod(commands , lib, origProt , 2, randomFactorSCMOD);
    }

	System.out.println("************* Minimizing the hydrogens of reference *********");
	distanceMatrix = new DistanceMatrix(copyReference.atoms(), 5.5, 2.0,4); 
	energy = new TotalEnergy(copyReference, distanceMatrix, energyCreatorsNoInf, commands);
	minimizer = new LBFGS(energy, commands);  
	minimizer.minimize();
	copyReference.defrost();
	
	
	System.out.println("******************** Phase II - free minimization ********************");
	System.out.println("************* Energy of the model *********");
	/* Initial RMS,GDT */
   	protein = new Protein(origProt.atoms().duplicate() , new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
  	protein.defrost();
	System.out.println("OUT: 111111 00 " + GDTcalculator.gdt(reference.atoms(),protein.atoms()));	 
	System.out.println("OUT: 222222 00 " + reference.atoms().CAFilter().getRms(protein.atoms().CAFilter()));
	/* Initial energies of the model */ 
	distanceMatrix = new DistanceMatrix(protein.atoms(), 5.5, 2.0, 4);  
	energy = new TotalEnergy(protein, distanceMatrix, energyCreatorsNoInf, commands);
	energy.evaluate();
	SolvateEnergy solvTerm = (SolvateEnergy) energy.getEnergyTerm(new SolvateEnergy());
	CompositePropensity2DEnergy propTerm = (CompositePropensity2DEnergy) energy.getEnergyTerm(new CompositePropensity2DEnergy());
	System.out.println("OUT: 333333 00 " + energy.report(2) + " " + 
						solvTerm.evaluate(false,1.0,1.0,1.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,0.0,1.0) + " 0.0");
	System.out.println("OUT: 999999 "+iteration);
	// PDB OUT    protein.atoms().print(new MeshiWriter(refFileName+".00." + modelNum + ".decoy"));
	// Same thing now with tether 
	System.out.println("************* Minimizing the model with short thether *********");
   	protein = new Protein(origProt.atoms().duplicate() , new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
  	protein.defrost();
	distanceMatrix = new DistanceMatrix(protein.atoms(), 5.5, 2.0, 4);  
	energy = new TotalEnergy(protein, distanceMatrix, energyCreatorsNoInf, commands);
	minimizer = new LBFGS(energy, 0.05, 300, 100);  
	System.out.println(minimizer.minimize());
	System.out.println("OUT: 777777 0 " + GDTcalculator.gdt(reference.atoms(),protein.atoms()));	 
	System.out.println("OUT: 888888 0 " + reference.atoms().CAFilter().getRms(protein.atoms().CAFilter()));
	solvTerm = (SolvateEnergy) energy.getEnergyTerm(new SolvateEnergy());
	propTerm = (CompositePropensity2DEnergy) energy.getEnergyTerm(new CompositePropensity2DEnergy());
	System.out.println("OUT: 888777 0 " + energy.report(2) + " " + 
						solvTerm.evaluate(false,1.0,1.0,1.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,0.0,1.0) + " " +
						origProt.atoms().CAFilter().getRms(protein.atoms().CAFilter()));
	System.out.println("OUT: 999999 "+iteration); 
   	origProt = new Protein(protein.atoms().duplicate() , new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
	// PDB OUT   protein.atoms().print(new MeshiWriter(refFileName+".0." + modelNum + ".decoy"));

    tmpEnergy = Wprop*propTerm.evaluate() + Wsolv*solvTerm.evaluate(false,1.0,1.0,1.0,0.0) + Whb*solvTerm.evaluate(false,0.0,0.0,0.0,1.0);
    if (tmpEnergy<bestEnergy) {
    	bestEnergy = tmpEnergy;
    	bestList = iteration;
    	chosenList = protein.atoms().duplicate();
    }
    do {
    	System.out.println("************ Inflating **********"); 
    	protein = new Protein(origProt.atoms().duplicate() , new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
    	protein.defrost();
    	protein.atoms().renumber();
		distanceMatrix = new DistanceMatrix(protein.atoms(), 5.5, 2.0, 4);  
		energy = new TotalEnergy(protein, distanceMatrix, energyCreators, commands);
		minimizer = new LBFGS(energy, 0.05, 200000, 100);  
		System.out.println(minimizer.minimize());
    	System.out.println("************ Minimizing **********");     	
		distanceMatrix = new DistanceMatrix(protein.atoms(), 5.5, 2.0, 4);  
		energy = new TotalEnergy(protein, distanceMatrix, energyCreatorsNoInf, commands);
		TetherEnergy tetherTerm = (TetherEnergy) energy.getEnergyTerm(new TetherEnergy());
		tetherTerm.off();
		solvTerm = (SolvateEnergy) energy.getEnergyTerm(new SolvateEnergy());
		propTerm = (CompositePropensity2DEnergy) energy.getEnergyTerm(new CompositePropensity2DEnergy());
		minimizer = new LBFGS(energy, 0.05, 20000, 100);  
		System.out.println(minimizer.minimize());
	   	iteration++;
		System.out.println("OUT: 111111 " + iteration + " " + GDTcalculator.gdt(reference.atoms(),protein.atoms()));	 
		System.out.println("OUT: 222222 " + iteration + " " + reference.atoms().CAFilter().getRms(protein.atoms().CAFilter()));
		System.out.println("OUT: 333333 " + iteration + " " + energy.report(2) + " " + 
							solvTerm.evaluate(false,1.0,1.0,1.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,0.0,1.0) + " " +  
							origProt.atoms().CAFilter().getRms(protein.atoms().CAFilter()));
		System.out.println("OUT: 999999 "+iteration); 
    	tmpEnergy = Wprop*propTerm.evaluate() + Wsolv*solvTerm.evaluate(false,1.0,1.0,1.0,0.0) + Whb*solvTerm.evaluate(false,0.0,0.0,0.0,1.0);
    	if (tmpEnergy<bestEnergy) {
    		bestEnergy = tmpEnergy;
	    	bestList = iteration;
    		chosenList = protein.atoms().duplicate();
	    }
//	   	if (iteration<numOfSCMODiterations)
//	   		SCMOD.scmod(commands , lib, protein , 2, randomFactorSCMOD);  
	// PDB OUT   protein.atoms().print(new MeshiWriter(refFileName+"." + iteration + "." + modelNum + ".decoy"));
	} while (iteration<numOfSCMODiterations);


	if (modelNum==1) {
	System.out.println("******************** Phase III - minimization of the reference ********************");
	 /* Initial RMS,GDT */
	distanceMatrix = new DistanceMatrix(copyReference.atoms(), 5.5, 2.0, 4);  
	energy = new TotalEnergy(copyReference, distanceMatrix, energyCreatorsNoInf, commands);
	solvTerm = (SolvateEnergy) energy.getEnergyTerm(new SolvateEnergy());
	energy.evaluate();
	System.out.println("OUT: 444444 0 " + GDTcalculator.gdt(reference.atoms(),copyReference.atoms()));	 
	System.out.println("OUT: 555555 0 " + reference.atoms().CAFilter().getRms(copyReference.atoms().CAFilter()));
	System.out.println("OUT: 666666 0 " + energy.report(2) + " " + 
			solvTerm.evaluate(false,1.0,1.0,1.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,0.0,1.0));
	minimizer = new LBFGS(energy, 0.05, 300, 100);  
	System.out.println(minimizer.minimize());
	
	System.out.println("OUT: 444444 0 " + GDTcalculator.gdt(reference.atoms(),copyReference.atoms()));	 
	System.out.println("OUT: 555555 0 " + reference.atoms().CAFilter().getRms(copyReference.atoms().CAFilter()));
	System.out.println("OUT: 666666 0 " + energy.report(2) + " " + 
			solvTerm.evaluate(false,1.0,1.0,1.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,0.0,1.0));
	TetherEnergy tetherTerm = (TetherEnergy) energy.getEnergyTerm(new TetherEnergy());
	tetherTerm.off();
	minimizer = new LBFGS(energy, 0.05, 20000, 100);  
	System.out.println(minimizer.minimize());

	System.out.println("OUT: 444444 1 " + GDTcalculator.gdt(reference.atoms(),copyReference.atoms()));	 
	System.out.println("OUT: 555555 1 " + reference.atoms().CAFilter().getRms(copyReference.atoms().CAFilter()));
	System.out.println("OUT: 666666 1 " + energy.report(2) + " " + 
			solvTerm.evaluate(false,1.0,1.0,1.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,0.0,1.0));
	// PDB OUT   copyReference.atoms().print(new MeshiWriter(refFileName+".000." + modelNum + ".decoy"));
	}

	/** OUTPUT - The final minimized coordinates are printed to the standard output. Note that 
	 * protein.atoms() returns the atom list of the protein we minimized. The print() method 
	 * of the list prints then atoms in a standard pdb format.
	 **/
	 System.out.println("OUT: 000000 " + bestList);
	 chosenList.print();  
    }



	protected static AtomList getMatchingAtoms(AtomList refAtoms , AtomList protAtoms) {
		AtomList result = new AtomList();
		for (int c=0 ; c<refAtoms.size() ; c++) 
		try {
			result.add(protAtoms.findAtomInList(refAtoms.atomAt(c).name , refAtoms.atomAt(c).residueNumber()));
		}
		catch (Exception e) {
			System.out.println("\n\n" + refAtoms.atomAt(c) + "\n\n");
			throw new RuntimeException(e);
		}
		return result;
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


	String line;
	String errorMessage = ("\n                  ******************\n"+
			       "Usage java -Xmx300m MinimizeProtein <commands file name> <pdb file name> seed\n"+
			       "                    ******************\n");
			      
	if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));
	commandsFileName = getOrderedArgument(args);
	if (commandsFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# commandsFileName = "+commandsFileName);

	commands = new CommandList(commandsFileName);
	
	modelFileName = getOrderedArgument(args);
	if (modelFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# initial model file name is "+modelFileName);

	refFileName = getOrderedArgument(args);
	if (refFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# reference file name is "+refFileName);


	String seedString = getOrderedArgument(args);
	if (seedString== null) throw new RuntimeException(errorMessage);
	int seed = (new Integer(seedString)).intValue();
	System.out.println("# seed is "+seed);
	initRandom(seed);


	String tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	modelNum = (new Integer(tmpString)).intValue();
	System.out.println("# modelNum is " + modelNum);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	numOfSCMODiterations = (new Integer(tmpString)).intValue();
	System.out.println("# numOfSCMODiterations is "+numOfSCMODiterations);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	withSCMOD = (new Integer(tmpString)).intValue();
	System.out.println("# withSCMOD (true if >0) is "+withSCMOD);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	randomFactorSCMOD = (new Double(tmpString)).doubleValue();
	System.out.println("# randomFactorSCMOD is " + randomFactorSCMOD);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	Wrg = (new Double(tmpString)).doubleValue();
	System.out.println("# Wrg is " + Wrg);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	Wev = (new Double(tmpString)).doubleValue();
	System.out.println("# Wev is " + Wev);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	Wsolv = (new Double(tmpString)).doubleValue();
	System.out.println("# Wsolv is " + Wsolv);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	Whb = (new Double(tmpString)).doubleValue();
	System.out.println("# Whb is " + Whb);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	Wprop = (new Double(tmpString)).doubleValue();
	System.out.println("# Wprop is " + Wprop);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	Wramach = (new Double(tmpString)).doubleValue();
	System.out.println("# Wramach is " + Wramach);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	Wsc = (new Double(tmpString)).doubleValue();
	System.out.println("# Wsc is " + Wsc);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	Winf = (new Double(tmpString)).doubleValue();
	System.out.println("# Winf is " + Winf);

	tmpString = getOrderedArgument(args);
	if (tmpString== null) throw new RuntimeException(errorMessage);
	infTar = (new Double(tmpString)).doubleValue();
	System.out.println("# infTar is " + infTar);
    }
}


