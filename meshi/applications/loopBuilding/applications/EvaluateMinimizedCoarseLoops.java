package meshi.applications.loopBuilding.applications;

import java.text.DecimalFormat;
import java.util.StringTokenizer;

import programs.PutHydrogens;

import meshi.applications.corpus.PropensityMatrix;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.ROT1solvation.CentroidSolvationCreator;
import meshi.energy.ROT1solvation.ROT1SolvationEnergy;
import meshi.energy.compositeTorsions.ramachandran.RamachandranCreator;
import meshi.energy.compositeTorsions.ramachandran.RamachandranEnergy;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBondEnergy;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_HighAccuracy_BBonly_Creator;
import meshi.energy.softExcludedVol.SoftExcludedVol;
import meshi.energy.softExcludedVol.SoftExcludedVolCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.ResidueBuilder;
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
 * This program will evaluate the energies of minimized loops. The usage is:
 *
 * java -Xmx300m EvaluateMinimizedCoarseLoops <commands file name> <PDB file name of Ref> <PDB file name of template> <file with list of loop models> <file with list of phipsi data> 
 * <loop starting resisue> <loop ending residue> <0 - for cluster evaluation ; 1 - for singlton evaluation>
 * 
 * The output is printed to the screen. 'grep' the lines that start with '999999'. The format of one of these lines is: 
 *
 * 999999 <0-based number in the list> <RMS-BB-global> <RMS-all-heavy-global> <total energy> <energy term 1> <energy term 2> .... 
 * 
 **/

public class EvaluateMinimizedCoarseLoops extends MeshiProgram implements Residues, AtomTypes{ 

	private static CommandList commands; 
	private static String commandsFileName = null;
	private static String modelFileName = null;  
	private static String loopsFileName = null;  
	private static String dataFileName = null;  
	private static String refFileName = null;  
	private static boolean isSingle = false;  
	private static int resStart = -999;
	private static int resEnd = -999;


	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		Protein reference = null;
		Protein model = null;
		DistanceMatrix distanceMatrix = null;
		TotalEnergy energy = null;
		DecimalFormat fmt2 = new DecimalFormat("0.###");

		// for 'noClust' energy selection on 8mers
		double Wcent_8 = 1.00;
		double Wrama_8 = 0.5;
		double Whb_8 = 1.0;
		double Wev_8 = 4.0;
		double Wprop_8 = 0.1;
		double Wclu_8 = -0.0;

		// for 'noClust' energy selection on 12mers                                                                                                                                                              
		double Wcent_12 = 1.0;
		double Wrama_12 = 0.0;
		double Whb_12 = 1.25;
		double Wev_12 = 3.0;
		double Wprop_12 = 0.4;
		double Wclu_12 = -0.0;

		// for Cluster energy selection on 8mers
		double Wcent_c_8 = 1.00;
		double Wrama_c_8 = 0.0;
		double Whb_c_8 = 0.9;
		double Wev_c_8 = 3.3;
		double Wprop_c_8 = 0.45;
		double Wclu_c_8 = -0.5;

		// for Cluster energy selection on 12mers                                                                                                                                                              
		double Wcent_c_12 = 1.0;
		double Wrama_c_12 = 0.0;
		double Whb_c_12 = 0.5;
		double Wev_c_12 = 1.85;
		double Wprop_c_12 = 0.1;
		double Wclu_c_12 = -0.5;

		// for 'noClust' energy selection	
		double Wcent = Double.NaN;
		double Wrama = Double.NaN;
		double Whb = Double.NaN;
		double Wev = Double.NaN;
		double Wprop = Double.NaN;
		double Wclu = Double.NaN;

		if (isSingle) {
			if ((resEnd-resStart+1)>=12) {
				Wcent = Wcent_12;
				Wrama = Wrama_12;
				Whb = Whb_12;
				Wev = Wev_12;
				Wprop = Wprop_12;
				Wclu = Wclu_12*(resEnd-resStart+1)/12.0;
			}
			else if ((resEnd-resStart+1)>8) {
				Wcent = Wcent_8 + (Wcent_12-Wcent_8)*((resEnd-resStart+1)-8)/4.0;
				Wrama = Wrama_8 + (Wrama_12-Wrama_8)*((resEnd-resStart+1)-8)/4.0;
				Whb = Whb_8 + (Whb_12-Whb_8)*((resEnd-resStart+1)-8)/4.0;
				Wev = Wev_8 + (Wev_12-Wev_8)*((resEnd-resStart+1)-8)/4.0;
				Wprop = Wprop_8 + (Wprop_12-Wprop_8)*((resEnd-resStart+1)-8)/4.0;
				Wclu = Wclu_8 + (Wclu_12-Wclu_8)*((resEnd-resStart+1)-8)/4.0;
			}
			else if ((resEnd-resStart+1)>3) {
				Wcent = Wcent_8;
				Wrama = Wrama_8;
				Whb = Whb_8;
				Wev = Wev_8;
				Wprop = Wprop_8;
				Wclu = Wclu_8*(resEnd-resStart+1)/8.0;			
			} 
			else
				throw new RuntimeException("Currently handling only more than 3 mers loops.");
		}
		else {
			if ((resEnd-resStart+1)>=12) {
				Wcent = Wcent_c_12;
				Wrama = Wrama_c_12;
				Whb = Whb_c_12;
				Wev = Wev_c_12;
				Wprop = Wprop_c_12;
				Wclu = Wclu_c_12*(resEnd-resStart+1)/12.0;
			}
			else if ((resEnd-resStart+1)>8) {
				Wcent = Wcent_c_8 + (Wcent_c_12-Wcent_c_8)*((resEnd-resStart+1)-8)/4.0;
				Wrama = Wrama_c_8 + (Wrama_c_12-Wrama_c_8)*((resEnd-resStart+1)-8)/4.0;
				Whb = Whb_c_8 + (Whb_c_12-Whb_c_8)*((resEnd-resStart+1)-8)/4.0;
				Wev = Wev_c_8 + (Wev_c_12-Wev_c_8)*((resEnd-resStart+1)-8)/4.0;
				Wprop = Wprop_c_8 + (Wprop_c_12-Wprop_c_8)*((resEnd-resStart+1)-8)/4.0;
				Wclu = Wclu_c_8 + (Wclu_c_12-Wclu_c_8)*((resEnd-resStart+1)-8)/4.0;
			}
			else if ((resEnd-resStart+1)>3) {
				Wcent = Wcent_c_8;
				Wrama = Wrama_c_8;
				Whb = Whb_c_8;
				Wev = Wev_c_8;
				Wprop = Wprop_c_8;
				Wclu = Wclu_c_8*(resEnd-resStart+1)/8.0;			
			} 
			else
				throw new RuntimeException("Currently handling only more than 3 mers loops.");	
		}

		
		EnergyCreator[] energyCreators = {  
//				new BondCreator(),
//				new AngleCreator(),
//				new PlaneCreator(),
//				new OutOfPlaneCreator(),
				new SoftExcludedVolCreator(1.0 , 12 , 1.0),
				new SimpleHydrogenBond_Dahiyat_HighAccuracy_BBonly_Creator(1.0),
				new RamachandranCreator(1.0),
				new CentroidSolvationCreator(1.0,"8.0",false)
		};	


		// Loading the reference and the background model	
//		reference = new Protein((new AtomList(refFileName)).noOXTFilter().backbone(), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
//		model = new Protein((new AtomList(modelFileName)).noOXTFilter().backbone(), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		reference = new Protein((new AtomList(refFileName)).noOXTFilter().backbone(), new ResidueExtendedAtoms(ADD_HYDROGENS_AND_FREEZE));
		PutHydrogens.adjustHydrogens(commands, reference);
		model = new Protein((new AtomList(modelFileName)).noOXTFilter().backbone(), new ResidueExtendedAtoms(ADD_HYDROGENS_AND_FREEZE));
		PutHydrogens.adjustHydrogens(commands, model);
		reference.defrost();
		model.defrost();
		for (int res=0 ; res<reference.residues().size(); res++) {
			if ((reference.residues().residueAt(res).type<20) &&
					(reference.residues().residueAt(res).type>-1)) {
				try {
					ResidueBuilder.buildCentroid(reference.residues().residueAt(res));
				}
				catch (Exception ee) {}
			}
		}		
		for (int cc=0 ; cc<model.atoms().size() ; cc++)
			model.atoms().atomAt(cc).setChain("A");
		double[][] originalPositions = new double[model.atoms().size()][3];
		for (int cc=0 ; cc<model.atoms().size() ; cc++) {
			originalPositions[cc][0] = model.atoms().atomAt(cc).x();
			originalPositions[cc][1] = model.atoms().atomAt(cc).y();
			originalPositions[cc][2] = model.atoms().atomAt(cc).z();
		}
		
		
		
		distanceMatrix = new DistanceMatrix(model.atoms(), 8.0, 1.0, 4);  
		energy = new TotalEnergy(model, distanceMatrix, energyCreators, commands);
		SimpleHydrogenBondEnergy hbEnergy = (SimpleHydrogenBondEnergy) energy.getEnergyTerm(new SimpleHydrogenBondEnergy());
		RamachandranEnergy ramachEnergy = (RamachandranEnergy) energy.getEnergyTerm(new RamachandranEnergy());
		ROT1SolvationEnergy CentTerm = (ROT1SolvationEnergy) energy.getEnergyTerm(new ROT1SolvationEnergy());
		SoftExcludedVol evTerm = (SoftExcludedVol) energy.getEnergyTerm(new SoftExcludedVol());
		PropensityMatrix propMatrix = new PropensityMatrix(commands,"allResPropensity.dat");
				

		// Looping on all the models
		String[] models = File2StringArray.f2a(loopsFileName);
		String[] dataFile = File2StringArray.f2a(dataFileName);

		for (int i=0 ; i<models.length ; i++) {
			try {	
				AtomList loop = (new AtomList(models[i])).backbone();
				for (int c=0; c<loop.size() ; c++) {
					Atom atom = loop.atomAt(c);
					model.atoms().findAtomInList(atom.name(), atom.residueNumber()).setXYZ(atom.x(), atom.y(), atom.z());
				}

				// Setting the minimization energy
				String[] phipsiData = File2StringArray.f2a(dataFile[i]);
				StringTokenizer st = new StringTokenizer(phipsiData[0]);
				st.nextToken();
				st.nextToken();
				st.nextToken();
				double clusterSize = Double.parseDouble(st.nextToken());

				energy.update();
				double[][] pp = RotamericTools.phipsi(model, distanceMatrix);
				double Ehb = hbEnergy.evaluate();
				double Eramach = ramachEnergy.evaluate();
				double Eev = evTerm.evaluate();
				double Eprop = calcProp(pp,propMatrix);
				for (int res=0 ; res<model.residues().size(); res++) {
					if ((model.residues().residueAt(res).type<20) &&
							(model.residues().residueAt(res).type>-1)) {
						try {
							ResidueBuilder.buildCentroid(model.residues().residueAt(res));
						}
						catch (Exception ee) {}
					}
				}		
				energy.update();
				double Ecent = CentTerm.evaluate();
				double rmsAll = calcRMSallHeavyAtoms(model, reference, resStart, resEnd);
				double rmsBB = calcRMS(model, reference, resStart, resEnd);
				for (int cc=0 ; cc<model.atoms().size() ; cc++) {
					model.atoms().atomAt(cc).setXYZ(originalPositions[cc][0], 
							originalPositions[cc][1], 
							originalPositions[cc][2]);
				}
				double totalEnergyVal = Wcent*Ecent + Wrama*Eramach + Whb*Ehb + Wev*Eev + Wprop*Eprop + Wclu*Math.log(1.0+clusterSize);				

				System.out.println("999999 " + i + " " + fmt2.format(rmsBB) + " " + fmt2.format(rmsAll) + " " + fmt2.format(totalEnergyVal) + " " +
						fmt2.format(Ecent) + " " + fmt2.format(Eramach) + " " + fmt2.format(Ehb) + " " +
						fmt2.format(Eev) + " " + fmt2.format(Eprop) + " " + clusterSize/*fmt2.format(Math.log(clusterSize+1))*/);
			}
			catch (Exception e) {
				System.out.println("Evaluation was not successful.\n");
				e.printStackTrace();
				System.out.println("\nContinueing to next model.\n");
			}
		}

	} // Of main


	
	private static double calcProp(double[][] pp, PropensityMatrix propMatrix) {
		double sum = 0.0;
		for (int res=resStart-1 ; res<=(resEnd+1) ; res++) {
			sum += propMatrix.propVal((int) pp[res][2], pp[res][0], pp[res][1]);			
		}
		return sum;
	}
	

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
			atom = prot.residue(c).atoms().findAtomInList("O",c);
			atomr = ref.residue(c).atoms().findAtomInList("O",c);
			totRms += (atom.x() - atomr.x())*
			(atom.x() - atomr.x()) +
			(atom.y() - atomr.y())*
			(atom.y() - atomr.y()) + 
			(atom.z() - atomr.z())*
			(atom.z() - atomr.z());
		}
		return Math.sqrt(totRms/(4*(end-start+1)));
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


		String errorMessage = ("\n                  ******************\n"+
				"Usage java -Xmx300m EvaluateMinimizedCoarseLoops <commands file name> <PDB file name of Ref> <PDB file name of template> <file with list of loop models> <file with list of phipsi data> " +
				"<loop starting resisue> <loop ending residue> <0 - for cluster evaluation ; 1 - for singlton evaluation>\n"+
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));
		commandsFileName = getOrderedArgument(args);
		if (commandsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# commandsFileName = "+commandsFileName);

		commands = new CommandList(commandsFileName);

		refFileName = getOrderedArgument(args);
		if (refFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# reference file name is "+refFileName);

		modelFileName = getOrderedArgument(args);
		if (modelFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# initial model file name is "+modelFileName);

		loopsFileName = getOrderedArgument(args);
		if (loopsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Loops loops file name is "+loopsFileName);

		dataFileName = getOrderedArgument(args);
		if (dataFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# data file name is "+dataFileName);

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
		if (Integer.parseInt(tmpString)==0) {
			isSingle = false;
			System.out.println("# This is a cluster evaluation.");
		}
		else if (Integer.parseInt(tmpString)==1) {
			isSingle = true;
			System.out.println("# This is a singleton evaluation.");
		}
		else {
			throw new RuntimeException(errorMessage);
		}
	}
}




