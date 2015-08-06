package programs;

import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.CoulombElectrostatics.CoulombElectrostaticCreator;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.compositeTorsions.compositePropensity2DwithPP.CompositePropensity2DCreator;
import meshi.energy.compositeTorsions.ramachandranAndChi1.RamachandranAndChi1PartialResiduesCreator;
import meshi.energy.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.softExcludedVol.SoftExcludedVolCreator;
import meshi.energy.solvate.SolvateCreatorRegularHB;
import meshi.energy.solvate.SolvateEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;

/**
 *<pre>
 * This program will minimize a batch of loops given in a list, and will compare the results to a reference.
 * Only the loop residues are compared in the RMS. There is no superposition, and it is assumed the rest of
 * the protein is aligned with the reference.
 * Unix usage:
 *     java -Xmx300m MinimizeBatchOfLoops <commands file name> <file with list of pdbs> <ref pdb file name> <loop starting resisue> <loop ending residue> <Wrg> <Wev> <Wsolv> <Whb> <Wprop> <Wramach>
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

public class EvaluateBatchOfLoops extends MeshiProgram implements Residues, AtomTypes{ 

    private static CommandList commands; 
    private static String commandsFileName = null;
    private static String modelsFileName = null;  
    private static String refFileName = null;  
    private static int resStart = -999;
    private static int resEnd = -999;
    private static double Wrg = 0.0;  
    private static double Wev = 3.0;  
    private static double Wsolv = 0.5;  
    private static double Whb = 1.0;  
    private static double Wprop = 1.0;  
    private static double Wramach = 0.1;  

 
    public static void main(String[] args) {
	init(args); 
	Protein reference = null;
	Protein model = null;
	DistanceMatrix distanceMatrix = null;
	TotalEnergy energy = null;
	SolvateEnergy solvTerm = null;

	// The creators for the terms
	EnergyCreator[] energyCreators = {  
		    new BondCreator(),
		    new AngleCreator(),
		    new PlaneCreator(),
		    new OutOfPlaneCreator(),
		    new CoulombElectrostaticCreator(1.0),
 /* I changed this line on 3/1/2008 to an electrostatic creator:
 			new LinearRgCreator(Wrg),   */
			new SoftExcludedVolCreator(Wev , 3 , 1.0),
//			new SolvateCreatorHBforMinimization(Wsolv,Wsolv,Wsolv,Whb),
			new SolvateCreatorRegularHB(Wsolv,Wsolv,Wsolv,Whb),
			new CompositePropensity2DCreator(Wprop),
			new RamachandranSidechainEnergyCreator(Wramach),
			new RamachandranAndChi1PartialResiduesCreator(Wramach)
		};	
	
	
//	// The creators for the terms
//	EnergyCreator[] energyCreators = {  
//		    new BondCreator(),
//		    new AngleCreator(),
//		    new PlaneCreator(),
//		    new OutOfPlaneCreator(),
//		    new CoulombElectrostaticCreator(1.0),
//			new EVenergyCreator(Wev , 0 , 1.0),
//			new SolvateCreatorRegularHB(0.001,Wsolv,0.001,Whb),
//			new CompositePropensity2DCreator(Wprop),
//			new RamachandranSidechainEnergyCreator(Wramach),
//			new RamachandranAndChi1PartialResiduesCreator(Wramach),
//		};	
	
	// Loading the reference	
	reference = new Protein(new AtomList(refFileName), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
	
	// Model files
	String[] models = File2StringArray.f2a(modelsFileName);

	// Reading the model for the first time
	model = new Protein(new AtomList(models[0]) , new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
	for (int cc=0 ; cc<model.atoms().size() ; cc++)
			model.atoms().atomAt(cc).setChain("A");
	
	// Energy and RMS 
    model.freeze();
    for (int c=resStart ; c<=resEnd ; c++)
    	model.residue(c).atoms().defrost();
	distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
	energy = new TotalEnergy(model, distanceMatrix, energyCreators, commands);
	energy.evaluate();
	solvTerm = (SolvateEnergy) energy.getEnergyTerm(new SolvateEnergy());
	
	// Looping on all the models	
	for (int i=0 ; i<models.length ; i++) {
		Protein tmpProt = new Protein(new AtomList(models[i].replace(".min.pdb", "")) , new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
	    for (int c=resStart ; c<=resEnd ; c++)
	    	for (int d=0 ; d<model.residue(c).atoms().size() ; d++) {
	    		Atom atom = tmpProt.residue(c).atoms().findAtomInList(model.residue(c).atoms().atomAt(d).name(), c);
	    		model.residue(c).atoms().atomAt(d).setXYZ(atom.x(), atom.y(), atom.z());    		
	    	}
		energy.evaluate();
		System.out.println("999999 " + i + " 0000 " + models[i].replace(".min.pdb", ""));
		System.out.println("999999 " + i + " 1111 " + calcRMS(model, reference, resStart, resEnd) + " " +
				calcRMSallHeavyAtoms(model, reference, resStart, resEnd));	 
		/* Initial energies of the model */ 
		System.out.println("999999 " + i + " 2222 " + energy.report(2) + " " + 
							solvTerm.evaluate(false,1.0,0.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,1.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,1.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,0.0,1.0));
		// Energy and RMS - After minimization
		tmpProt = new Protein(new AtomList(models[i]) , new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
	    for (int c=resStart ; c<=resEnd ; c++)
	    	for (int d=0 ; d<model.residue(c).atoms().size() ; d++) {
	    		Atom atom = tmpProt.residue(c).atoms().findAtomInList(model.residue(c).atoms().atomAt(d).name(), c);
	    		model.residue(c).atoms().atomAt(d).setXYZ(atom.x(), atom.y(), atom.z());    		
	    	}
		energy.evaluate();
		System.out.println("999999 " + i + " 3333 " + calcRMS(model, reference, resStart, resEnd) + " " +
				calcRMSallHeavyAtoms(model, reference, resStart, resEnd));	 
		/* Final energies of the model */ 
		System.out.println("999999 " + i + " 4444 " + energy.report(2) + " " + 
				solvTerm.evaluate(false,1.0,0.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,1.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,1.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,0.0,1.0));
		energy.evaluateAtoms();
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
    			       "Usage java -Xmx300m MinimizeBatchOfLoops <commands file name> <file with list of pdbs> <ref pdb file name> " +
    			       "<loop starting resisue> <loop ending residue> <Wrg>  <Wev> <Wsolv> <Whb> <Wprop> <Wramach> \n"+
    			       "                    ******************\n");
    			      
    	if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));
    	commandsFileName = getOrderedArgument(args);
    	if (commandsFileName == null) throw new RuntimeException(errorMessage);
    	System.out.println("# commandsFileName = "+commandsFileName);

    	commands = new CommandList(commandsFileName);
    	
    	modelsFileName = getOrderedArgument(args);
    	if (modelsFileName == null) throw new RuntimeException(errorMessage);
    	System.out.println("# initial model file name is "+modelsFileName);

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
        }
}
