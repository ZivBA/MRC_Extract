import alignment.RvalAlignerCluster;
import utils.ExtractMaxValue;
import utils.ScoreUtilities.ScoringGeneralHelpers;
import utils.Scoring.MRC_Score;
import utils.UtilExceptions.MissingChainID;

import java.io.File;
import java.io.IOException;

import java.util.Arrays;

/**
 * worker thread class - each instance processes an input chain from start to finish.
 * creating a worker instance queues the chains on the executor service to be run when resources are available
 */
public class WorkerThread implements Runnable {

	private String[] args;

	public WorkerThread(String[] s) {
		this.args = s;
	}

	@Override
	public void run() {

		MRC_Score myScore = null;
		try {
			// if argument string does not include chain to process, then entire protein needs to be processed (not used).
			if (args.length == 2) {
				myScore = new MRC_Score(args[1], args[0]);
				char[] chains = myScore.getMyProt().getChains();


				for (char chain : chains){
					myScore = new MRC_Score(args[1], args[0], String.valueOf(chain),myScore.getMyMap());
					myScore.requestedChain = chain;
					myScore.scoreProtein();
					myScore.calcZvalue();
					myScore.createCSVs();
					System.out.println("Done with: " + args[0] + " chain: " + myScore.requestedChain);
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
					String swissProtPath = "/home/zivben/IdeaProjects/Results/uniprot_sprot_10_2015.fasta";

					String seqListPath = myScore.getMyProt().getSource().getParent() +File.separator + filePrefix + ".fasta";
					//			String profileFilePath1 = folderPath + filePrefix + "_profileNoVec.txt";
					String profileFilePath2 = folderPath + filePrefix + "_profileLatestVec.txt";
					String profileFilePath3 = folderPath + filePrefix + "_profileNoVec_weighteBB.txt";
					String profileWeight2 = folderPath + filePrefix + "_profileLatestVec_weightedBB_2.txt";
					String profileWeight5 = folderPath + filePrefix + "_profileLatestVec_weightedBB_5.txt";
					String profileWeight10 = folderPath + filePrefix + "_profileLatestVec_weightedBB_10.txt";

					//			RvalAlignerCluster.runThread(swissProtPath, seqListPath, profileFilePath1);
//					RvalAlignerCluster.runThread(swissProtPath, seqListPath, profileFilePath2);
//					//			RvalAlignerCluster.runThread(swissProtPath, seqListPath, profileFilePath3);
//					RvalAlignerCluster.runThread(swissProtPath, seqListPath, profileWeight2);
//					RvalAlignerCluster.runThread(swissProtPath, seqListPath, profileWeight5);
//					RvalAlignerCluster.runThread(swissProtPath, seqListPath, profileWeight10);
				}


			// if arg list includes chain to process, ignore other chains and process just that one.
			} else if (args.length == 3) {
				myScore = new MRC_Score(args[1], args[0], args[2]);     // Create MRC map object. in the process a SimpleProtein object is created.
				float[] maxValResult = ExtractMaxValue.getMaxValue(myScore.getMyMap());
				System.out.println(Arrays.toString(maxValResult));
				ExtractMaxValue.writeMarkerFile(myScore.getMyProt().getSource().getParent(), maxValResult);
				myScore.scoreProtein();
				myScore.calcZvalue();
				myScore.createCSVs();
				System.out.println("Done with: " + args[0] + " Chain " + args[2]);
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
				String filePrefix = myScore.getMyProt().getFileName() + "_" + myScore.requestedChain;
				String swissProtPath = "/home/zivben/IdeaProjects/Results/uniprot_sprot_10_2015.fasta";

				String seqListPath = myScore.getMyProt().getSource().getParent() +File.separator + filePrefix + ".fasta";
				//			String profileFilePath1 = folderPath + filePrefix + "_profileNoVec.txt";
				String profileFilePath2 = folderPath + filePrefix + "_profileLatestVec.txt";
				String profileFilePath3 = folderPath + filePrefix + "_profileNoVec_weighteBB.txt";
				String profileWeight2 = folderPath + filePrefix + "_profileLatestVec_weightedBB_2.txt";
				String profileWeight5 = folderPath + filePrefix + "_profileLatestVec_weightedBB_5.txt";
				String profileWeight10 = folderPath + filePrefix + "_profileLatestVec_weightedBB_10.txt";

				//			RvalAlignerCluster.runThread(swissProtPath, seqListPath, profileFilePath1);
				//			RvalAlignerCluster.runThread(swissProtPath, seqListPath, profileFilePath2);
				//			RvalAlignerCluster.runThread(swissProtPath, seqListPath, profileFilePath3);
				RvalAlignerCluster.runThread(swissProtPath, seqListPath, profileWeight2);
				RvalAlignerCluster.runThread(swissProtPath, seqListPath, profileWeight5);
				RvalAlignerCluster.runThread(swissProtPath, seqListPath, profileWeight10);
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
