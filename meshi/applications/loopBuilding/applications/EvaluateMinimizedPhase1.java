package meshi.applications.loopBuilding.applications;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Vector;

import meshi.applications.loopBuilding.AbstractLoopBuilder;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.ROT1solvation.CBSolvationCreator;
import meshi.energy.compositeTorsions.compositePropensity2DwithPP.CompositePropensity2DCreator;
import meshi.energy.compositeTorsions.ramachandran.RamachandranCreator;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_HighAccuracy_BBonly_Creator;
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
import meshi.util.clustering.Cluster;
import meshi.util.clustering.HierarchicalClusterer;
import meshi.util.file.File2StringArray;
import meshi.util.file.MeshiWriter;
import meshi.util.rotamericTools.RotamericTools;



public class EvaluateMinimizedPhase1 extends MeshiProgram implements Residues, AtomTypes{ 

	private static CommandList commands; 
	private static String commandsFileName = null;
	private static String modelsFileName = null;  
	private static String modelFileName = null;  
	private static String refFileName = null;  
	private static int resStart = -999;
	private static int resEnd = -999;
	private static int Noutputs = -999;
	private static String outputPath = null;  
	private static final double rmsClusteringCutoff = 1.5;
	private static final int upperSizeClusteringCutoff = 50;
	private static final int lowerSizeClusteringCutoff = 2;


	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		DunbrackLib lib = new DunbrackLib(commands, 1.0, 1);
		Protein reference = null;
		Protein model = null;
		DistanceMatrix distanceMatrix = null;
		TotalEnergy energy = null;
		DecimalFormat fmt2 = new DecimalFormat("0.###");

		// Loading the proteins for clustering
		String[] models = File2StringArray.f2a(modelsFileName);
		Protein[] loops = new Protein[models.length];
		for (int i=0 ; i<models.length ; i++) {
			loops[i] = new Protein((new AtomList(models[i])).noOXTFilter(), new ResidueExtendedAtoms(ADD_ATOMS));
			for (int cc=0 ; cc<loops[i].atoms().size() ; cc++)
				loops[i].atoms().atomAt(cc).setChain("A");
			RotamericTools.putIntoRot1(loops[i], new DistanceMatrix(loops[i].atoms(), 5.5, 2.0, 3), lib);
		}
		double[][] disMat = new double[models.length][models.length];
		for (int i=0 ; i<models.length ; i++) {
			for (int j=i+1 ; j<models.length ; j++) {
				disMat[i][j] = calcRMSallHeavyAtoms(loops[i], loops[j], resStart, resEnd);
				disMat[j][i] = disMat[i][j];
			}
		}
		
		// Clustering
		HierarchicalClusterer clusterer = new HierarchicalClusterer(disMat);
		clusterer.cluster(rmsClusteringCutoff, upperSizeClusteringCutoff);

		// Calculating the energies
		EnergyCreator[] energyCreatorsFinalEvaluation = {  
				new SimpleHydrogenBond_Dahiyat_HighAccuracy_BBonly_Creator(1.0),
				new RamachandranCreator(0.2),
				new CompositePropensity2DCreator(0.6),
				new CBSolvationCreator(1.5,"9.0",false,null)
		};
		reference = new Protein((new AtomList(refFileName)).noOXTFilter(), new ResidueExtendedAtoms(ADD_ATOMS));
		for (int cc=0 ; cc<reference.atoms().size() ; cc++)
			reference.atoms().atomAt(cc).setChain("A");
		RotamericTools.putIntoRot1(reference, new DistanceMatrix(reference.atoms(), 5.5, 2.0, 3), lib);
		model = new Protein((new AtomList(modelFileName)).noOXTFilter(), new ResidueExtendedAtoms(ADD_ATOMS));
		for (int cc=0 ; cc<model.atoms().size() ; cc++)
			model.atoms().atomAt(cc).setChain("A");
		RotamericTools.putIntoRot1(model, new DistanceMatrix(model.atoms(), 5.5, 2.0, 3), lib);
		model.freeze();
		for (int cc=0 ; cc<loops[0].atoms().size() ; cc++)
			if (loops[0].atoms().atomAt(cc).name().equals("CA"))
				loops[0].atoms().atomAt(cc).residue().atoms().defrost();
		distanceMatrix = new DistanceMatrix(model.atoms(), 9.0, 0.1, 3);  
		energy = new TotalEnergy(model, distanceMatrix, energyCreatorsFinalEvaluation, commands);
		double[] eneVal = new double[models.length];
		double[] rmsBB = new double[models.length];
		double[] rmsAllAtoms = new double[models.length];
		for (int i=0 ; i<models.length ; i++) {
			// Reading the model
			AtomList loop = new AtomList(models[i]);
			for (int c=0; c<loop.size() ; c++) {
				Atom atom = loop.atomAt(c);
				model.atoms().findAtomInList(atom.name(), atom.residueNumber()).setXYZ(atom.x(), atom.y(), atom.z());
			}
			energy.evaluate(); // Just so that the distance matrix is updated
			RotamericTools.putIntoRot1(model, RotamericTools.phipsi(model, distanceMatrix), lib);
			eneVal[i] = energy.evaluate();
			rmsBB[i] = calcRMS(model, reference, resStart, resEnd);
			rmsAllAtoms[i] = calcRMSallHeavyAtoms(model, reference, resStart, resEnd); 
		}
		
		// analyzing the clusters with the energies
		Vector <Cluster> clusters = new Vector <Cluster>();
		Vector <Double> clustEne = new Vector <Double>();
		Vector <Integer> clustCenter = new Vector <Integer>();
		clusterer.initializeSerialTokenizer();
		for (Cluster clust=clusterer.getNextSerialToken() ; (clust!=null) && (clust.getClusterMembers().size()>1) ; clust=clusterer.getNextSerialToken()) {
			if (clust.getClusterMembers().size()>lowerSizeClusteringCutoff) {
				clusters.add(clust);
				double sumEne = 0.0;
				double sumNum = 0.0;
				double minDis = Double.MAX_VALUE;
				double sumDis = 0.0;
				int minInd = -1;
				for (int cc=0 ; cc<clust.getClusterMembers().size() ; cc++) {
					sumEne += eneVal[clust.getClusterMembers().get(cc).intValue()];
					sumNum++;
					sumDis = 0.0;
					for (int cc1=0 ; cc1<clust.getClusterMembers().size() ; cc1++) {
						sumDis += disMat[clust.getClusterMembers().get(cc).intValue()][clust.getClusterMembers().get(cc1).intValue()];
					}
					if (sumDis<minDis) {
						minDis = sumDis;
						minInd = clust.getClusterMembers().get(cc).intValue();
					}
				}
				clustEne.add(new Double(sumEne/sumNum - 0.05*clust.getClusterMembers().size()));
				clustCenter.add(new Integer(minInd));
			}
		}
		
		if (clusters.size()==0) {
			System.out.println("NoClustError: No meaningful clusters were found.");			
		}
		else {
			for (int cc=0 ; cc<clusters.size() ; cc++) {
				System.out.println("Cluster: " + cc + "   with Energy: " + clustEne.get(cc).doubleValue() + " and center at: " + clustCenter.get(cc).intValue());
				for (int ccc=0; ccc<clusters.get(cc).getClusterMembers().size() ; ccc++) {
					int ind = clusters.get(cc).getClusterMembers().get(ccc).intValue();
					System.out.println(ind + " " + models[ind] + " " + rmsBB[ind] + " " + rmsAllAtoms[ind] + " " + eneVal[ind]);
				}
			}
		}

		
		// Outputing:
		double[] clusterEnergies = new double[clusters.size()];
		for (int cc=0 ; cc<clusters.size() ; cc++) {
			clusterEnergies[cc] = clustEne.get(cc).doubleValue();
		}
		int[] sortedClusterIndices = AbstractLoopBuilder.findTopMinArray(clusterEnergies, clusters.size(), Double.MAX_VALUE);
		for (int cc=0 ; (cc<clusters.size()) && (cc<Noutputs) ; cc++) {
			int ind = clustCenter.get(sortedClusterIndices[cc]).intValue();
			System.out.println("Printing to disk model " + cc + " that came from:" + models[ind] + " of cluster:" + sortedClusterIndices[cc] + " and the following RMS (BB and all-atom):" + rmsBB[ind] + " " + rmsAllAtoms[ind]);
			try {
				loops[ind].atoms().noOXTFilter().print(new MeshiWriter(outputPath+"/"+cc+".pdb"));
			} catch (IOException e) {
				System.out.println("Could not write the output: " + outputPath+"/"+cc+".pdb");
				e.printStackTrace();
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
				if (!prot.residue(c).atoms().atomAt(d).isHydrogen && !prot.residue(c).atoms().atomAt(d).equals("OXT")) {
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
				"Usage java -Xmx300m MinimizeBatchOfCoarseLoops <commands file name> <file with list of pdbs> <file with list of phipsi data> <ref pdb file name> " +
				"<loop starting resisue> <loop ending residue> <Wev> <Whb> <Wtorval> <Wtether> <output PDB extension string>\n"+
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

		modelsFileName = getOrderedArgument(args);
		if (modelsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Loop models file name is "+modelsFileName);

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
		Noutputs = (new Integer(tmpString)).intValue();
		System.out.println("# Number of output models: " + Noutputs);
		
		outputPath = getOrderedArgument(args);
		if (outputPath == null) throw new RuntimeException(errorMessage);
		System.out.println("# The output path: " + outputPath);		
	}
}
