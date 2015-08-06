package programs;

import java.util.StringTokenizer;

import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.compositeTorsions.compositePropensity2DwithPP.CompositePropensity2DCreator;
import meshi.energy.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.excludedVolumeImprovedDistance.EVenergyCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.solvateNew.SolvateCreatorRegularHB;
import meshi.energy.solvateNew.SolvateEnergy;
import meshi.energy.tether.TetherCreator;
import meshi.energy.tether.TetherEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
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

/**
 *<pre>
 *
 * Unix usage:
 *     java -Xmx300m MinimizeLoopGridSearch <commands file name> <model file name> <ref pdb file name> <loop starting resisue> <loop ending residue> <Weights file> 
 *
 **/

public class MinimizeLoopGridSearch extends MeshiProgram implements Residues, AtomTypes{ 

    private static CommandList commands; 
    private static String commandsFileName = null;
    private static String weightsFileName = null;  
    private static String modelFileName = null;  
    private static String refFileName = null;  
    private static String endString = null;
    private static int resStart = -999;
    private static int resEnd = -999;
    private static double Wrg = 0.0;  
    private static double Wev = 3.0;  
    private static double Wsolv = 0.5;  
    private static double Whb = 1.0;  
    private static double Wprop = 1.0;  
    private static double Wramach = 0.1;  

 
    public static void main(String[] args) throws MinimizerException, LineSearchException{
	init(args); 
	Protein reference = null;
	Protein model = null;
	DistanceMatrix distanceMatrix = null;
	TotalEnergy energy = null;
	SolvateEnergy solvTerm = null;
	Minimizer minimizer = null;

	// The creators for the terms
	EnergyCreator[] energyCreators = {  
		    new BondCreator(),
		    new AngleCreator(),
		    new PlaneCreator(),
		    new OutOfPlaneCreator(),
		    new EVenergyCreator(0.0,0,1.0),
//			new SoftExcludedVolCreator(Wev , 3 , 1.0),
		    new SolvateCreatorRegularHB(0.0001,0.0001,0.0001,0.0),
			new CompositePropensity2DCreator(0.0),
			new RamachandranSidechainEnergyCreator(0.0),
			new TetherCreator(5.0, new AtomList.ClassCaCbFilter())
		};	
	
	// Loading the reference	
	reference = new Protein(new AtomList(refFileName), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));

	// Loading the model	
	model = new Protein(new AtomList(modelFileName), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
	double[][] initialcoordinates = saveLoopCoordinates(model);  
	for (int cc=0 ; cc<model.atoms().size() ; cc++)
		model.atoms().atomAt(cc).setChain("A");
	model.freeze();
	for (int c=resStart ; c<=resEnd ; c++)
		model.residue(c).atoms().defrost();

	// Looping on all the models
	String[] weightStrings = File2StringArray.f2a(weightsFileName);
	for (int i=0 ; i<weightStrings.length ; i++) {
		restoreLoopCoordinates(model, initialcoordinates);
		double[] weights = str2vec(weightStrings[i]);
		String header = "999999" + i + "9999 ";
		energyCreators[4].setWeight(weights[0]); // EV
		((SolvateCreatorRegularHB) energyCreators[5]).setHBweight(weights[1]); // HB
		((SolvateCreatorRegularHB) energyCreators[5]).setCarbonWeight(weights[2]); // Carbon
		energyCreators[6].setWeight(weights[3]); // Prop
		energyCreators[7].setWeight(weights[4]); // Ramach
	  	
		System.out.println(header + i + " 0000 " + i);
		System.out.println(header + i + " 1111 " + calcRMS(model, reference, resStart, resEnd) + " " +
				calcRMSallHeavyAtoms(model, reference, resStart, resEnd));	 
		/* Initial energies of the model */ 
		distanceMatrix = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);  
		energy = new TotalEnergy(model, distanceMatrix, energyCreators, commands);
		energy.evaluate();
		solvTerm = (SolvateEnergy) energy.getEnergyTerm(new SolvateEnergy());
		System.out.println(header + i + " 2222 " + energy.report(2) + " " + 
							solvTerm.evaluate(false,1.0,0.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,1.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,1.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,0.0,1.0));

		// Energy and RMS - After minimization
		minimizer = new LBFGS(energy, 0.05, 750, 500);
		try {
			System.out.println(minimizer.minimize());
		}
		catch (Exception e) {
			System.out.println("ERROR in minimization");
		}
		(energy.getEnergyTerm(new TetherEnergy())).off();
		System.out.println("Tether is OFF");
		minimizer = new LBFGS(energy, 0.05, 20000, 500);
		try {
			System.out.println(minimizer.minimize());
		}
		catch (Exception e) {
			System.out.println("ERROR in minimization");
		}
		System.out.println(header + i + " 3333 " + calcRMS(model, reference, resStart, resEnd) + " " +
				calcRMSallHeavyAtoms(model, reference, resStart, resEnd));	 
		/* Final energies of the model */ 
		energy.evaluate();
		System.out.println(header + i + " 4444 " + energy.report(2) + " " + 
				solvTerm.evaluate(false,1.0,0.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,1.0,0.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,1.0,0.0) + " " + solvTerm.evaluate(false,0.0,0.0,0.0,1.0));
		energy.evaluateAtoms();
//		try {
//			model.atoms().print(new MeshiWriter(models[i]+ "." + endString));
//		}
//		catch (Exception e) {
//			System.out.print("\nThere was a problem writing the output:\n" + e + "\n\nContinueing...\n\n");
//		}
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

    private static double[] str2vec(String s) {
    	StringTokenizer st = new StringTokenizer(s);
    	double[] result = new double[st.countTokens()];
    	for (int c=0; st.countTokens()>0 ; c++)
    		result[c]=(new Double(st.nextToken())).doubleValue();
    	return result;
    }
    
    private static double[][] saveLoopCoordinates(Protein model) {
    	double[][] savedLoopCoordinates = new double[3][(resEnd-resStart+1)*15];
    	Atom atom;
    	int counter=0;
    	for (int c=resStart ; c<=resEnd ; c++) 
    		for (int cc=0 ; cc<model.residue(c).atoms().size() ; cc++){
    			atom = model.residue(c).atoms().atomAt(cc);
    			savedLoopCoordinates[0][counter] = atom.x();			
    			savedLoopCoordinates[1][counter] = atom.y();			
    			savedLoopCoordinates[2][counter] = atom.z();
    			counter++;
    		}	
    	return savedLoopCoordinates;
    }

    private static void restoreLoopCoordinates(Protein model, double[][] savedLoopCoordinates) {
    	Atom atom;
    	int counter=0;
    	for (int c=resStart ; c<=resEnd ; c++) 
    		for (int cc=0 ; cc<model.residue(c).atoms().size() ; cc++){
    			atom = model.residue(c).atoms().atomAt(cc);
    			atom.setXYZ(savedLoopCoordinates[0][counter],savedLoopCoordinates[1][counter],
    			savedLoopCoordinates[2][counter]);
    			counter++;
    		}	
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
    			       "Usage java -Xmx300m MinimizeLoopGridSearch <commands file name> <model file name> <ref pdb file name> <loop starting resisue> <loop ending residue> <Weights file>\n"+
    			       "                    ******************\n");
    			      
    	if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));
    	commandsFileName = getOrderedArgument(args);
    	if (commandsFileName == null) throw new RuntimeException(errorMessage);
    	System.out.println("# commandsFileName = "+commandsFileName);

    	commands = new CommandList(commandsFileName);
    	
    	modelFileName = getOrderedArgument(args);
    	if (modelFileName == null) throw new RuntimeException(errorMessage);
    	System.out.println("# loop model file name is "+modelFileName);

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
    	
    	weightsFileName = getOrderedArgument(args);
    	if (weightsFileName == null) throw new RuntimeException(errorMessage);
    	System.out.println("# Weights file name is "+weightsFileName);
    	
        }
}
