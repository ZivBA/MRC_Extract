package meshi.applications.loopBuilding;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.DecimalFormat;
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

public class LoopBuilderDevelopePhase2SecondAttempt extends
LoopBuilderUserDefinedFragmentsStochasticWithClosure {
	
	// User-defined CONST:
	// *******************
	// for demi-clustering
	protected final double COVERAGE_CRITERIA = 0.5; 
	protected final int CLUSTERS_FOR_COVERAGE_CRITERIA = 5; 
	protected final double CUTOFF_INCREAMENTS = 0.5;
	// for final clustering
	protected final double FRACTION_OF_LOOPS_FOR_MINIMAL_CLUSTER = 0.001; 
	protected final int NUMBER_OF_CLUSTERS_IN_OUTPUT = 50; 
	protected final int NUMBER_OF_SINGLE_IN_OUTPUT = 30; 

	// for 'noClust' energy selection on 8mers
	protected final double Wcent_8 = 1.00;
	protected final double Whbsr_8 = 0.6;
	protected final double Whblr_8 = 2.35;
	protected final double Wev_8 = 2.0;
	protected final double Wprop_8 = 0.65;
	protected final double Wclu_8 = -0.55;

	// for 'noClust' energy selection on 12mers                                                                                                                                                              
	protected final double Wcent_12 = 1.0;
	protected final double Whbsr_12 = 0.45;
	protected final double Whblr_12 = 0.9;
	protected final double Wev_12 = 0.25;
	protected final double Wprop_12 = 0.1;
	protected final double Wclu_12 = -0.15;

	protected double Wcent = Double.NaN;
	protected double Whbsr = Double.NaN;
	protected double Whblr = Double.NaN;
	protected double Wev = Double.NaN;
	protected double Wprop = Double.NaN;
	protected double Wclu = Double.NaN;
    
    
    
	protected ROT1SolvationEnergy cent_080 = null;
	protected SimpleHydrogenBondEnergy hb_sr = null;
	protected SimpleHydrogenBondEnergy hb_lr = null;
	
	protected double[][] intial_HM_LoopCoordinates = null;

	public LoopBuilderDevelopePhase2SecondAttempt(CommandList commands,
			String writePath, Corpus corpus, Protein prot, Protein ref,
			int resStart, int resEnd, double rmsMatchCO, double rmsCutOff,
			Vector<int[]> fragsDescription, double closureTolerance) {
		super(commands, writePath, corpus, prot, ref, resStart, resEnd,
				rmsMatchCO, rmsCutOff, fragsDescription, closureTolerance);
		if ((resEnd-resStart+1)>9) {
			Wcent = Wcent_12;
			Whbsr = Whbsr_12;
			Whblr = Whblr_12;
			Wev = Wev_12;
			Wprop = Wprop_12;
			Wclu = Wclu_12*(resEnd-resStart+1)/12.0;
		}
		else if ((resEnd-resStart+1)>3) {
			Wcent = Wcent_8;
			Whbsr = Whbsr_8;
			Whblr = Whblr_8;
			Wev = Wev_8;
			Wprop = Wprop_8;
			Wclu = Wclu_8*(resEnd-resStart+1)/8.0;			
		} 
		else
			throw new RuntimeException("Currently handling only more than 3 mers loops.");
		
		intial_HM_LoopCoordinates = saveLoopCoordinates();
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

/*		// Printing the distance matrix for adelene
		try{
			DecimalFormat fmt = new DecimalFormat("0.##");
			BufferedWriter bw = new BufferedWriter(new FileWriter(writePath+"/distanceMatrix.txt"));
			for (int i=0 ; i<allResults.size() ; i++) {
				for (int j=i+1 ; j<allResults.size() ; j++) {
					bw.write(i + " " + j + " " + fmt.format(disMat[i][j]) + "\n");
				}
			}
			bw.close();
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}		
*/	
		
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
		// Printing what was the initial RMS
		restoreLoopCoordinates(intial_HM_LoopCoordinates);
		System.out.println("\nThe initial BB-RMSD is: " + calcRMSonHonigBackbone(resStart, resEnd));
		
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
			array_Eev[modelNum] = allResults.get(modelNum).evEnergy;
			array_Eprop[modelNum] = allResults.get(modelNum).propEnergy;			
		}
		// ***************************
		// Building the energy vectors - End
		// ***************************
		
		
		// clustering and outputing
		double clustCutoff = 3.0;//findClusteringCriteria(disMat , COVERAGE_CRITERIA);
		clusteringAndOutputing(array_Ecent, array_Ehb_sr,
				array_Ehb_lr, array_Eev,
				array_Eprop, array_RMS_honig, 
				array_RMS_allHeavy, disMat, 
				clustCutoff, "777755", "666655");		
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
				double average_Ecent = clust.findPercentileOfTrait(array_Ecent,0.25);
				double average_Ehb_sr = clust.findPercentileOfTrait(array_Ehb_sr,0.25);
				double average_Ehb_lr = clust.findPercentileOfTrait(array_Ehb_lr,0.25);
				double average_Eev = clust.findPercentileOfTrait(array_Eev,0.25);
				double average_Eprop = clust.findPercentileOfTrait(array_Eprop,0.25);
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
			if (clusterAffiliation[modelNum]==-1) {
				System.out.print(headerAllLoops + " " + modelNum + " " +
						fmt2.format(array_RMS_honig[modelNum]) + " " +
						fmt2.format(array_RMS_allHeavy[modelNum]) + " " +
						fmt2.format(array_Ecent[modelNum]) + " " +
						fmt2.format(array_Ehb_sr[modelNum]) + " " +
						fmt2.format(array_Ehb_lr[modelNum]) + " " +
						fmt2.format(array_Eev[modelNum]) + " " +
						fmt2.format(array_Eprop[modelNum]) + "   " +
						fmt2.format(Math.log(allResults.size()*FRACTION_OF_LOOPS_FOR_MINIMAL_CLUSTER)) + " " +
						fmt2.format(allResults.size()*FRACTION_OF_LOOPS_FOR_MINIMAL_CLUSTER) + " " + 
						fmt2.format(clusterAffiliation[modelNum]) + "          ");		
			}
			else {
				System.out.print(headerAllLoops + " " + modelNum + " " +
						fmt2.format(array_RMS_honig[modelNum]) + " " +
						fmt2.format(array_RMS_allHeavy[modelNum]) + " " +
						fmt2.format(array_Ecent[modelNum]) + " " +
						fmt2.format(array_Ehb_sr[modelNum]) + " " +
						fmt2.format(array_Ehb_lr[modelNum]) + " " +
						fmt2.format(array_Eev[modelNum]) + " " +
						fmt2.format(array_Eprop[modelNum]) + "   " +
						fmt2.format(Math.log(clusterData.get(clusterAffiliation[modelNum])[5])) + " " +
						fmt2.format(clusterData.get(clusterAffiliation[modelNum])[5]) + " " + 
						fmt2.format(clusterAffiliation[modelNum]) + "          ");				
			}
			fragRank = allResults.get(modelNum).fragRank;
			for (int tmpc=0 ; tmpc<fragRank.length ; tmpc++) {
				System.out.print(fragRank[tmpc]+" ");
				//score += 0.0*(Math.log(1+fragRank[tmpc]));
			}
			System.out.println();
		}
		
		// Writing the lower 'noClust' loops
		double[] singleScores = new double[allResults.size()];
		for (int modelNum=0 ; modelNum<allResults.size() ; modelNum++) {
			if (clusterAffiliation[modelNum]==-1) {
				singleScores[modelNum] = Wcent*array_Ecent[modelNum] + 
				Whbsr*array_Ehb_sr[modelNum] + 
				Whblr*array_Ehb_lr[modelNum] + 
				Wev*array_Eev[modelNum] + 
				Wprop*array_Eprop[modelNum]+
				Wclu*Math.log(allResults.size()*FRACTION_OF_LOOPS_FOR_MINIMAL_CLUSTER);
			}
			else {
				singleScores[modelNum] = Wcent*array_Ecent[modelNum] + 
				Whbsr*array_Ehb_sr[modelNum] + 
				Whblr*array_Ehb_lr[modelNum] + 
				Wev*array_Eev[modelNum] + 
				Wprop*array_Eprop[modelNum]+
				Wclu*Math.log(clusterData.get(clusterAffiliation[modelNum])[5]);				
			}
		}
		int[] sortedSingleScores = AbstractLoopBuilder.findTopMinArray(singleScores, NUMBER_OF_SINGLE_IN_OUTPUT, Double.MAX_VALUE);
		try{
			for (int writeLoop = 0; writeLoop<sortedSingleScores.length ; writeLoop++) {
				System.out.println("Writing the single loop to disk: The index is " + 
						sortedSingleScores[writeLoop] + " with energy " + singleScores[sortedSingleScores[writeLoop]]);
				restoreLoopCoordinates(allResults.get(sortedSingleScores[writeLoop]).coors);
				DecimalFormat fmt = new DecimalFormat("0.##");
				BufferedWriter bw = new BufferedWriter(new FileWriter(writePath+"/single/"+writeLoop+".pdb"));
				for (int res=resStart ; res<=resEnd ; res++) {
					for (int atInd=0 ; atInd<prot.residue(res).atoms().size() ; atInd++) {
						bw.write(prot.residue(res).atoms().atomAt(atInd).toString() + "\n");
					}
				}
				bw.close();
				bw = new BufferedWriter(new FileWriter(writePath+"/single/"+writeLoop+".data"));
				for (int res=resStart ; res<=resEnd ; res++) {
					bw.write(res + " " + 
							fmt.format(allResults.get(sortedSingleScores[writeLoop]).myPhiPsi[res-resStart][0]) + " " + 
							fmt.format(allResults.get(sortedSingleScores[writeLoop]).myPhiPsi[res-resStart][1]) + " 0 \n");
				}
				bw.close();
			}
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
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
					fmt2.format(clusterData.get(affiliationNumber)[6]) + " " +
					+affiliationNumber);
		}
		System.out.println("\n\n\n\n");

		// Writing the '55Clust' loops
		try{
			for (int cc=0 ; (cc<sortedClusterSize.length) && (cc<NUMBER_OF_CLUSTERS_IN_OUTPUT) ; cc++) {
				int affiliationNumber = sortedClusterSize[cc];			
				int center = (int) Math.round(clusterData.get(affiliationNumber)[6]);
				restoreLoopCoordinates(allResults.get(center).coors);	
				BufferedWriter bw = new BufferedWriter(new FileWriter(writePath+"/clust/"+ cc +".pdb"));
				for (int res=resStart ; res<=resEnd ; res++) {
					for (int atInd=0 ; atInd<prot.residue(res).atoms().size() ; atInd++) {
						bw.write(prot.residue(res).atoms().atomAt(atInd).toString() + "\n");
					}
				}
				bw.close();
				bw = new BufferedWriter(new FileWriter(writePath+"/clust/"+ cc +".data"));
				for (int res=resStart ; res<=resEnd ; res++) {
					bw.write(res + " " + 
							fmt2.format(allResults.get(center).myPhiPsi[res-resStart][0]) + " " + 
							fmt2.format(allResults.get(center).myPhiPsi[res-resStart][1]) + " " + 
							fmt2.format(clusterData.get(affiliationNumber)[5]) + "\n");
				}
				bw.close();			
			}
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}

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
