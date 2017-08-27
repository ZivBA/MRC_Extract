import alignment.RvalAlignerCluster;
import utils.ExtractMaxValue;
import utils.ScoreUtilities.ScoringGeneralHelpers;
import utils.Scoring.MRC_Score;
import utils.UtilExceptions.MissingChainID;
import utils.scwrlIntegration.SCWRLactions;

import java.io.File;
import java.io.IOException;

import java.util.Arrays;

/**
 * worker thread class - each instance processes an input chain from start to finish.
 * creating a worker instance queues the chains on the executor service to be run when resources are available
 */
public class WorkerThread implements Runnable {

	private String[] args;
	public static String swissProtPath = "/home/zivben/IdeaProjects/Results/uniprot_sprot_10_2015.fasta";
	public WorkerThread(String[] s) {
		this.args = s;
	}

	@Override
	public void run() {

		MRC_Score myScore = null;
		try {
			// if argument string does not include chain to process, then entire protein needs to be processed (not used).
//			System.out.println("Run SCWRL once on src protein to fix any inconsistencies with file");
			
			if (args.length == 2) {
//				System.out.println("Create MRC_Score object");
				myScore = new MRC_Score(args[1], args[0]);
				char[] chains = myScore.getMyProt().getChains();
				myScore.getMyProt().saveOriginalPositions();

				for (char chain : chains){
					myScore = new MRC_Score(myScore.getMyProt(), String.valueOf(chain),myScore.getMyMap());
//					System.out.println("Scoring Protein");
					myScore.scoreProtein();
//					System.out.println("Running calculations");
					myScore.calcZvalue();
//					System.out.println("Writing CSVs");
					myScore.createCSVs();
//					System.out.println("Done with: " + args[0] + " chain: " + myScore.requestedChain);
					if (ScoringGeneralHelpers.debug) {
						for (String line : myScore.logFile) {
							System.out.println(line);
						}
					} else {

						for (File path : myScore.toDelete) {
							delete(path);
						}
					}

					String folderPath = myScore.getMyProt().getSource().getParent() +File.separator+ "tempCSVs" +File.separator;
					String filePrefix = myScore.getMyProt().getFileName().toUpperCase() + "_" + myScore.requestedChain;
					

					String seqListPath = myScore.getMyProt().getSource().getParent() +File.separator + filePrefix + ".fasta";
					String profileNoVec0 = folderPath + filePrefix + "_profileNoVec_weightedBB_0.0.txt";
					String profileNoVec05 = folderPath + filePrefix + "_profileNoVec_weightedBB_0.5.txt";
					String profileNoVec1 = folderPath + filePrefix + "_profileNoVec_weightedBB_1.0.txt";
					String profileNoVec2 = folderPath + filePrefix + "_profileNoVec_weightedBB_2.0.txt";
					String profileNoVec4 = folderPath + filePrefix + "_profileNoVec_weightedBB_4.0.txt";
					String profileNoVec8 = folderPath + filePrefix + "_profileNoVec_weightedBB_8.0.txt";
					String profileNoVec12 = folderPath + filePrefix + "_profileNoVec_weightedBB_12.0.txt";
					String profileLatestVec0 = folderPath + filePrefix + "_profileLatestVec_weightedBB_0.0.txt";
					String profileLatestVec05 = folderPath + filePrefix + "_profileLatestVec_weightedBB_0.5.txt";
					String profileLatestVec1 = folderPath + filePrefix + "_profileLatestVec_weightedBB_1.0.txt";
					String profileLatestVec2 = folderPath + filePrefix + "_profileLatestVec_weightedBB_2.0.txt";
					String profileLatestVec4 = folderPath + filePrefix + "_profileLatestVec_weightedBB_4.0.txt";
					String profileLatestVec8 = folderPath + filePrefix + "_profileLatestVec_weightedBB_8.0.txt";
					String profileLatestVec12 = folderPath + filePrefix + "_profileLatestVec_weightedBB_12.0.txt";

					// first arg = null and not swissProtPath if we want to just score against ref FASTA
//					System.out.println("Running Alignment");
					RvalAlignerCluster.runThread(null, seqListPath, profileNoVec0);
					RvalAlignerCluster.runThread(null, seqListPath, profileNoVec05);
					RvalAlignerCluster.runThread(null, seqListPath, profileNoVec1);
					RvalAlignerCluster.runThread(null, seqListPath, profileNoVec2);
					RvalAlignerCluster.runThread(null, seqListPath, profileNoVec4);
					RvalAlignerCluster.runThread(null, seqListPath, profileNoVec8);
					RvalAlignerCluster.runThread(null, seqListPath, profileNoVec12);
					//					RvalAlignerCluster.runThread(swissProtPath, seqListPath, profileLatestVec);
					RvalAlignerCluster.runThread(null, seqListPath, profileLatestVec0);
					RvalAlignerCluster.runThread(null, seqListPath, profileLatestVec05);
					RvalAlignerCluster.runThread(null, seqListPath, profileLatestVec1);
					RvalAlignerCluster.runThread(null, seqListPath, profileLatestVec2);
					RvalAlignerCluster.runThread(null, seqListPath, profileLatestVec4);
					RvalAlignerCluster.runThread(null, seqListPath, profileLatestVec8);
					RvalAlignerCluster.runThread(null, seqListPath, profileLatestVec12);
				}


			// if arg list includes chain to process, ignore other chains and process just that one.
			} else if (args.length == 3) {
				myScore = new MRC_Score(args[1], args[0], args[2]);
				char[] chains = myScore.getMyProt().getChains();
				myScore.getMyProt().saveOriginalPositions();
			
				//					System.out.println("Scoring Protein");
				myScore.scoreProtein();
				//					System.out.println("Running calculations");
				myScore.calcZvalue();
				//					System.out.println("Writing CSVs");
				myScore.createCSVs();
				//					System.out.println("Done with: " + args[0] + " chain: " + myScore.requestedChain);
				if (ScoringGeneralHelpers.debug) {
					for (String line : myScore.logFile) {
						System.out.println(line);
					}
				} else {
					
					for (File path : myScore.toDelete) {
						delete(path);
					}
				}
				
				String folderPath = myScore.getMyProt().getSource().getParent() + File.separator + "tempCSVs" + File.separator;
				String filePrefix = myScore.getMyProt().getFileName().toUpperCase() + "_" + myScore.requestedChain;
				
				
				String seqListPath = myScore.getMyProt().getSource().getParent() + File.separator + filePrefix + ".fasta";
				String profileNoVec0 = folderPath + filePrefix + "_profileNoVec_weightedBB_0.0.txt";
				String profileNoVec05 = folderPath + filePrefix + "_profileNoVec_weightedBB_0.5.txt";
				String profileNoVec1 = folderPath + filePrefix + "_profileNoVec_weightedBB_1.0.txt";
				String profileNoVec2 = folderPath + filePrefix + "_profileNoVec_weightedBB_2.0.txt";
				String profileNoVec4 = folderPath + filePrefix + "_profileNoVec_weightedBB_4.0.txt";
				String profileNoVec8 = folderPath + filePrefix + "_profileNoVec_weightedBB_8.0.txt";
				String profileNoVec12 = folderPath + filePrefix + "_profileNoVec_weightedBB_12.0.txt";
				String profileLatestVec0 = folderPath + filePrefix + "_profileLatestVec_weightedBB_0.0.txt";
				String profileLatestVec05 = folderPath + filePrefix + "_profileLatestVec_weightedBB_0.5.txt";
				String profileLatestVec1 = folderPath + filePrefix + "_profileLatestVec_weightedBB_1.0.txt";
				String profileLatestVec2 = folderPath + filePrefix + "_profileLatestVec_weightedBB_2.0.txt";
				String profileLatestVec4 = folderPath + filePrefix + "_profileLatestVec_weightedBB_4.0.txt";
				String profileLatestVec8 = folderPath + filePrefix + "_profileLatestVec_weightedBB_8.0.txt";
				String profileLatestVec12 = folderPath + filePrefix + "_profileLatestVec_weightedBB_12.0.txt";
				
				// first arg = null and not swissProtPath if we want to just score against ref FASTA
				//					System.out.println("Running Alignment");
				RvalAlignerCluster.runThread(null, seqListPath, profileNoVec0);
				RvalAlignerCluster.runThread(null, seqListPath, profileNoVec05);
				RvalAlignerCluster.runThread(null, seqListPath, profileNoVec1);
				RvalAlignerCluster.runThread(null, seqListPath, profileNoVec2);
				RvalAlignerCluster.runThread(null, seqListPath, profileNoVec4);
				RvalAlignerCluster.runThread(null, seqListPath, profileNoVec8);
				RvalAlignerCluster.runThread(null, seqListPath, profileNoVec12);
				//					RvalAlignerCluster.runThread(swissProtPath, seqListPath, profileLatestVec);
				RvalAlignerCluster.runThread(null, seqListPath, profileLatestVec0);
				RvalAlignerCluster.runThread(null, seqListPath, profileLatestVec05);
				RvalAlignerCluster.runThread(null, seqListPath, profileLatestVec1);
				RvalAlignerCluster.runThread(null, seqListPath, profileLatestVec2);
				RvalAlignerCluster.runThread(null, seqListPath, profileLatestVec4);
				RvalAlignerCluster.runThread(null, seqListPath, profileLatestVec8);
				RvalAlignerCluster.runThread(null, seqListPath, profileLatestVec12);
				
				
			}


		} catch (IOException | MissingChainID e) {
			e.printStackTrace();
		}
	}

	private void delete(File path) {
		if (path.exists()) {
			File[] files = path.listFiles();
			for (int i = 0; i < files.length; i++) {
				if (files[i].isDirectory()) {
					delete(files[i]);
				} else {
					files[i].delete();
				}
			}
		}
		path.delete();
	}

}
