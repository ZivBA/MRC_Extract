package meshi.applications.loopBuilding;

import java.util.Vector;

import meshi.applications.corpus.Corpus;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.ROT1solvation.CentroidSolvationCreator;
import meshi.energy.ROT1solvation.ROT1SolvationEnergy;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBondEnergy;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_LowAccuracy_BBonly_LongRange_Creator;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_LowAccuracy_BBonly_ShortRange_Creator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.ResidueBuilder;
import meshi.molecularElements.Atom;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.clustering.Cluster;
import meshi.util.clustering.HierarchicalClusterer;

public class LoopBuilderDevelopePhase1SecondAttempt extends
LoopBuilderUserDefinedFragmentsStochasticWithClosure {
	
	// User-defined CONST:
	// *******************
	// for demi-clustering
	protected final int CLUSTERS_FOR_COVERAGE_CRITERIA = 5; 
	protected final double CUTOFF_INCREAMENTS = 0.5;
	// for final clustering
	protected final double FRACTION_OF_LOOPS_FOR_MINIMAL_CLUSTER = 0.001; 
	protected final int NUMBER_OF_CLUSTERS_IN_OUTPUT = 50; 
	
	protected ROT1SolvationEnergy cent_080 = null;
	protected SimpleHydrogenBondEnergy hb_sr = null;
	protected SimpleHydrogenBondEnergy hb_lr = null;

	public LoopBuilderDevelopePhase1SecondAttempt(CommandList commands,
			String writePath, Corpus corpus, Protein prot, Protein ref,
			int resStart, int resEnd, double rmsMatchCO, double rmsCutOff,
			Vector<int[]> fragsDescription, double closureTolerance) {
		super(commands, writePath, corpus, prot, ref, resStart, resEnd,
				rmsMatchCO, rmsCutOff, fragsDescription, closureTolerance);
	}

	protected void analyzedCloseLoop() {
		double evEnergy = evaluateEVextLoop()+evaluateEVinterLoop();
		double touch = evaluateTouchLoop();
		if ((evEnergy<EV_CUTOFF) && (touch>=((resEnd-resStart+1)/4))) {
			InternalLoopData results;
			try {
				results = evaluateScores(resStart ,resEnd); 	
			}
			catch (Exception e) {
				System.out.println("Warning: an error occured in the energy calculation. The loop will be dicarded.\n the error is:" + e.getMessage());
				return;
			}
			numberOfCoarseGenerated++;
			if (numberOfCoarseGenerated>maxNumberOfLoopGenerated)
				gotEnoughCoarseModels = true;
			if (numberOfCoarseGenerated/howOftenToPrintFinalLoop == ((int) numberOfCoarseGenerated/(1.0*howOftenToPrintFinalLoop))) {
				System.out.print(PRINTING_HEADER_COARSE + " " + numberOfCoarseGenerated + " " +
						fmt2.format(results.RMS) + " " +
						fmt2.format(evEnergy) + " " +
						fmt2.format(touch) + " " + 
						fmt2.format(results.bbHBenergy) + "         ");
				for (int tmpc=0 ; tmpc<fragRank.length ; tmpc++) {
					System.out.print(fragRank[tmpc]+" ");
					//score += 0.0*(Math.log(1+fragRank[tmpc]));
				}
				System.out.print("       ");
				for (int tmpc=0 ; tmpc<actuallyTaken.length ; tmpc++) {
					System.out.print(actuallyTaken[tmpc]+" ");
				}
				System.out.println();
			}
			allResults.add(new BasicLoopResult(evEnergy,results.propEnergy,results.bbHBenergy,0.0/*score is ignored*/,results.RMS,fragRank,saveLoopCoordinates(),saveLoopCoordinatesHeavyBackbone(),savePhiPsiOfLoop()));
		}
	}
	
	
	/** 
	 * WARNING: prot and ref are changed into their centroid rep!!!
	 */
	public void developementAnalysis(int moduleBase) {
		
		// Building distance matrix between loops
		double[][] disMat = new double[allResults.size()][allResults.size()];
		for (int i=0 ; i<allResults.size() ; i++) {
			for (int j=i+1 ; j<allResults.size() ; j++) {
				disMat[i][j] = allResults.get(i).calcRMS(allResults.get(j));
				disMat[j][i] = disMat[i][j];
			}
		}
		
		// ***************************
		// Building the energy vectors - Start
		// ***************************
		// The energy arrays 
		double[] array_Ecent = new double[allResults.size()];
		double[] array_Ehb_sr = new double[allResults.size()];
		double[] array_Ehb_lr = new double[allResults.size()];
		double[] array_Eev = new double[allResults.size()];
		double[] array_Eprop = new double[allResults.size()];
		double[] array_RMS_bb = new double[allResults.size()];
		double[] array_RMS_honig = new double[allResults.size()];
		double[] array_RMS_allHeavy = new double[allResults.size()];
		// Before doing the centroid representaion let's calculate Honig's RMS
		for (int modelNum=0 ; modelNum<allResults.size() ; modelNum++) {
			restoreLoopCoordinates(allResults.get(modelNum).coors);
			array_RMS_honig[modelNum] = calcRMSonHonigBackbone(resStart, resEnd);
		}
		// Putting the protein in centroid position
		for (int res=0 ; res<prot.residues().size(); res++) {
			if ((prot.residues().residueAt(res).type<20) &&
					(prot.residues().residueAt(res).type>-1))
				ResidueBuilder.buildCentroid(prot.residues().residueAt(res));
		}
		for (int res=0 ; res<ref.residues().size(); res++) {
			if ((ref.residues().residueAt(res).type<20) &&
					(ref.residues().residueAt(res).type>-1))
				ResidueBuilder.buildCentroid(ref.residues().residueAt(res));
		}
		// The creator and terms
		EnergyCreator[] energyCreators = {new CentroidSolvationCreator(1.0,"8.0",false),
				new SimpleHydrogenBond_Dahiyat_LowAccuracy_BBonly_LongRange_Creator(1.0,false),
				new SimpleHydrogenBond_Dahiyat_LowAccuracy_BBonly_ShortRange_Creator(1.0,false)
//				new SimpleHydrogenBond_Dahiyat_LowAccuracy_BBonly_InLoop_Creator(1.0,false,resStart,resEnd),
//				new SimpleHydrogenBond_Dahiyat_LowAccuracy_BBonly_Creator(1.0,false)
		};
		energy = new TotalEnergy(prot, new DistanceMatrix(prot.atoms(),  8.0, 0.1, 4), energyCreators, commands);
		cent_080 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[0]);
		hb_lr = (SimpleHydrogenBondEnergy) (energy.getEnergyTerms(new SimpleHydrogenBondEnergy())[0]);
		hb_sr = (SimpleHydrogenBondEnergy) (energy.getEnergyTerms(new SimpleHydrogenBondEnergy())[1]);
		// Looping on the loops
		for (int modelNum=0 ; modelNum<allResults.size() ; modelNum++) {
			restoreLoopCoordinates(allResults.get(modelNum).coors);
			for (int res=resStart ; res<=resEnd; res++) {
				ResidueBuilder.buildCentroid(prot.residue(res));
			}			
			try {
				energy.update();
			}
			catch (Exception e) {
				System.out.println("Warning: an error occured in the energy calculation of the development stage. The " + modelNum + " loop will be dicarded.\n the error is:" + e.getMessage());
			}
			array_RMS_bb[modelNum] = allResults.get(modelNum).rms;			
			array_RMS_allHeavy[modelNum] = calcRMSonHeavyBackbone(resStart, resEnd);
			array_Ecent[modelNum] = cent_080.evaluate();
			array_Ehb_sr[modelNum] = hb_sr.evaluate();
			array_Ehb_lr[modelNum] = hb_lr.evaluate();
//			array_Ehb_lr[modelNum] = hb_lr.evaluate() - array_Ehb_sr[modelNum];
			array_Eev[modelNum] = allResults.get(modelNum).evEnergy;
			array_Eprop[modelNum] = allResults.get(modelNum).propEnergy;			
		}
		// ***************************
		// Building the energy vectors - End
		// ***************************
		
		
		// clustering and outputing
//		clusteringAndOutputing(array_Ecent, array_Ehb_sr,
//				array_Ehb_lr, array_Eev,
//				array_Eprop, array_RMS_honig, 
//				array_RMS_allHeavy, disMat, 
//				1.0, "777710", "666610");		
//		clusteringAndOutputing(array_Ecent, array_Ehb_sr,
//				array_Ehb_lr, array_Eev,
//				array_Eprop, array_RMS_honig, 
//				array_RMS_allHeavy, disMat, 
//				1.5, "777715", "666615");		
		clusteringAndOutputing(array_Ecent, array_Ehb_sr,
				array_Ehb_lr, array_Eev,
				array_Eprop, array_RMS_honig, 
				array_RMS_allHeavy, disMat, 
				2.0, "777720", "666620");		
		clusteringAndOutputing(array_Ecent, array_Ehb_sr,
				array_Ehb_lr, array_Eev,
				array_Eprop, array_RMS_honig, 
				array_RMS_allHeavy, disMat, 
				2.5, "777725", "666625");		
//		clusteringAndOutputing(array_Ecent, array_Ehb_sr,
//				array_Ehb_lr, array_Eev,
//				array_Eprop, array_RMS_honig, 
//				array_RMS_allHeavy, disMat, 
//				3.0, "777730", "666630");		
//		clusteringAndOutputing(array_Ecent, array_Ehb_sr,
//				array_Ehb_lr, array_Eev,
//				array_Eprop, array_RMS_honig, 
//				array_RMS_allHeavy, disMat, 
//				4.0, "777740", "666640");		
		double clustCutoff = findClusteringCriteria(disMat , 0.5);
		clusteringAndOutputing(array_Ecent, array_Ehb_sr,
				array_Ehb_lr, array_Eev,
				array_Eprop, array_RMS_honig, 
				array_RMS_allHeavy, disMat, 
				clustCutoff, "777755", "666655");		
//		clustCutoff = findClusteringCriteria(disMat , 0.8);
//		clusteringAndOutputing(array_Ecent, array_Ehb_sr,
//				array_Ehb_lr, array_Eev,
//				array_Eprop, array_RMS_honig, 
//				array_RMS_allHeavy, disMat, 
//				clustCutoff, "777788", "666688");		
//		clustCutoff = findClusteringCriteria(disMat , 0.9);
//		clusteringAndOutputing(array_Ecent, array_Ehb_sr,
//				array_Ehb_lr, array_Eev,
//				array_Eprop, array_RMS_honig, 
//				array_RMS_allHeavy, disMat, 
//				clustCutoff, "777799", "666699");		
	}


	/** 
	 * WARNING: no actual energy calculations here. It's just for fragment picking assessment.
	 */
	public void fragmentPickingAnalysis(int moduleBase) {		
		// ***************************
		// Building the energy vectors - Start
		// ***************************
		// The energy arrays 
		double[] array_Eev = new double[allResults.size()];
		double[] array_Eprop = new double[allResults.size()];
		double[] array_RMS_bb = new double[allResults.size()];
		double[] array_RMS_honig = new double[allResults.size()];
		double[] array_RMS_allHeavy = new double[allResults.size()];
		// Before doing the centroid representaion let's calculate Honig's RMS
		for (int modelNum=0 ; modelNum<allResults.size() ; modelNum++) {
			restoreLoopCoordinates(allResults.get(modelNum).coors);
			array_RMS_honig[modelNum] = calcRMSonHonigBackbone(resStart, resEnd);
			restoreLoopCoordinates(allResults.get(modelNum).coors);
			array_RMS_bb[modelNum] = allResults.get(modelNum).rms;			
			array_RMS_allHeavy[modelNum] = calcRMSonHeavyBackbone(resStart, resEnd);
			array_Eev[modelNum] = allResults.get(modelNum).evEnergy;
			array_Eprop[modelNum] = allResults.get(modelNum).propEnergy;			
		}
		// ***************************
		// Building the energy vectors - End
		// ***************************

		// Outputing loop Data:
		for (int modelNum=0 ; modelNum<allResults.size() ; modelNum++) {
			// The data on the loop
			System.out.print("777777 " + modelNum + " " +
					fmt2.format(array_RMS_honig[modelNum]) + " " +
					fmt2.format(array_RMS_allHeavy[modelNum]) + " " +
					0.111 + " " +
					0.111 + " " +
					0.111 + " " +
					fmt2.format(array_Eev[modelNum]) + " " +
					fmt2.format(array_Eprop[modelNum]) + "   " +
					99 + "          ");
			fragRank = allResults.get(modelNum).fragRank;
			for (int tmpc=0 ; tmpc<fragRank.length ; tmpc++) {
				System.out.print(fragRank[tmpc]+" ");
				//score += 0.0*(Math.log(1+fragRank[tmpc]));
			}
			System.out.println();
		}
		
		System.out.println("\n\n\n\n");		
	}

	
	private double calcRMSonHonigBackbone(int start, int end) {
		double totRms = 0.0;
		int ntot = 0;
		for (int c=start; c<=end ; c++) {
			for (int d=0; d<prot.residue(c).atoms().size() ; d++) 
				if (!prot.residue(c).atoms().atomAt(d).isHydrogen &&
						prot.residue(c).atoms().atomAt(d).isBackbone &&
						!prot.residue(c).atoms().atomAt(d).name().equals("CB")) {
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
	
	
	private double calcRMSonHeavyBackbone(int start, int end) {
		double totRms = 0.0;
		int ntot = 0;
		for (int c=start; c<=end ; c++) {
			for (int d=0; d<prot.residue(c).atoms().size() ; d++) 
				if (!prot.residue(c).atoms().atomAt(d).isHydrogen &&
						prot.residue(c).atoms().atomAt(d).isBackbone) {
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
	
	public void setProt(Protein newProt) {
		prot = newProt;
	}
	public void setRef(Protein newRef) {
		ref = newRef;
	}


	
	private void clusteringAndOutputing(double[] array_Ecent, double[] array_Ehb_sr,
			double[] array_Ehb_lr, double[] array_Eev,
			double[] array_Eprop, double[] array_RMS_honig, 
			double[] array_RMS_allHeavy, double[][] disMat, 
			double rmsClusteringCutoff, String headerAllLoops, String headerClusters) {
		// ***************************
		// Clustering - Start
		// ***************************
		int[] clusterAffiliation = new int[allResults.size()];
		for (int modelNum=0 ; modelNum<allResults.size() ; modelNum++) {
			clusterAffiliation[modelNum] = -1;
		}
		HierarchicalClusterer clusterer = new HierarchicalClusterer(disMat);
		clusterer.cluster(rmsClusteringCutoff, Integer.MAX_VALUE);
		clusterer.initializeSerialTokenizer();
		Vector<double[]> clusterData = new Vector<double[]>(); 
		int clusterCounter = 0;
		for (Cluster clust=clusterer.getNextSerialToken() ; (clust!=null) && (clust.getClusterMembers().size()>1) ; clust=clusterer.getNextSerialToken()) {
			if (clust.getClusterMembers().size()>(allResults.size()*FRACTION_OF_LOOPS_FOR_MINIMAL_CLUSTER)) {
				double average_Ecent = clust.findMeanOfTrait(array_Ecent);
				double average_Ehb_sr = clust.findMeanOfTrait(array_Ehb_sr);
				double average_Ehb_lr = clust.findMeanOfTrait(array_Ehb_lr);
				double average_Eev = clust.findMeanOfTrait(array_Eev);
				double average_Eprop = clust.findMeanOfTrait(array_Eprop);
				for (int cc=0 ; cc<clust.getClusterMembers().size() ; cc++) {
					clusterAffiliation[clust.getClusterMembers().get(cc).intValue()] = clusterCounter;
				}
				double[] tmpArray = {average_Ecent, average_Ehb_sr, average_Ehb_lr, average_Eev, average_Eprop, 
						clust.getClusterMembers().size(), clust.findCenter()};
				clusterData.add(tmpArray);
				clusterCounter++;
			}
		}
		// ***************************
		// Clustering - End
		// ***************************

		// Outputing loop Data:
		for (int modelNum=0 ; modelNum<allResults.size() ; modelNum++) {
			// The data on the loop
			System.out.print(headerAllLoops + " " + modelNum + " " +
					fmt2.format(array_RMS_honig[modelNum]) + " " +
					fmt2.format(array_RMS_allHeavy[modelNum]) + " " +
					fmt2.format(array_Ecent[modelNum]) + " " +
					fmt2.format(array_Ehb_sr[modelNum]) + " " +
					fmt2.format(array_Ehb_lr[modelNum]) + " " +
					fmt2.format(array_Eev[modelNum]) + " " +
					fmt2.format(array_Eprop[modelNum]) + "   " +
					fmt2.format(clusterAffiliation[modelNum]) + "          ");
			fragRank = allResults.get(modelNum).fragRank;
			for (int tmpc=0 ; tmpc<fragRank.length ; tmpc++) {
				System.out.print(fragRank[tmpc]+" ");
				//score += 0.0*(Math.log(1+fragRank[tmpc]));
			}
			System.out.println();
		}
		
		// Outputting cluster Data:
		System.out.println("\nClusters:");
		double[] clusterSizes = new double[clusterData.size()];
		for (int cc=0 ; cc<clusterSizes.length ; cc++) {
			clusterSizes[cc] = -clusterData.get(cc)[5];
		}
		int[] sortedClusterSize = AbstractLoopBuilder.findTopMinArray(clusterSizes, clusterSizes.length, Double.MAX_VALUE);
		for (int cc=0 ; (cc<sortedClusterSize.length) && (cc<NUMBER_OF_CLUSTERS_IN_OUTPUT) ; cc++) {
			int affiliationNumber = sortedClusterSize[cc];			
			int center = (int) Math.round(clusterData.get(affiliationNumber)[6]);
			System.out.println(headerClusters + " " + cc + " " +
					fmt2.format(array_RMS_honig[center]) + " " +
					fmt2.format(array_RMS_allHeavy[center]) + " " +
					fmt2.format(clusterData.get(affiliationNumber)[0]) + " " +
					fmt2.format(clusterData.get(affiliationNumber)[1]) + " " +
					fmt2.format(clusterData.get(affiliationNumber)[2]) + " " +
					fmt2.format(clusterData.get(affiliationNumber)[3]) + " " +
					fmt2.format(clusterData.get(affiliationNumber)[4]) + "   " +
					fmt2.format(Math.log(clusterData.get(affiliationNumber)[5])) + " " +
					fmt2.format(clusterData.get(affiliationNumber)[5]) + "        " + 
					affiliationNumber);
		}
		System.out.println("\n\n\n\n");
	}	
	
	private double findClusteringCriteria(double[][] disMat , double coverage_criteria) {
		double initialRMScutOff = 0.0;
		double sumSize = 0;
		while (sumSize/allResults.size() < coverage_criteria) {
			initialRMScutOff += CUTOFF_INCREAMENTS;
			HierarchicalClusterer clusterer = new HierarchicalClusterer(disMat);
			clusterer.cluster(initialRMScutOff, Integer.MAX_VALUE);
			clusterer.initializeSizeTokenizer();
			sumSize = 0;
			int clusterCounter = 0;
			System.out.print("Summing clusters: ");
			for (Cluster clust=clusterer.getNextSizeToken() ; (clust!=null) && (clusterCounter<CLUSTERS_FOR_COVERAGE_CRITERIA) ; clusterCounter++) {
				System.out.print(clust.getClusterMembers().size() + " ");
				sumSize += clust.getClusterMembers().size();
				clust=clusterer.getNextSizeToken();
			}
			System.out.println();
		}
		System.out.println("The final RMS cutoff for clustering is: " + initialRMScutOff);		
		return initialRMScutOff;
	}
	
}
