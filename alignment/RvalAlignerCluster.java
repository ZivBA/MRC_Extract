package alignment;

import meshi.util.crossLinking.MySequence;
import meshi.util.crossLinking.MySequenceList;
import meshi.util.file.File2StringArray;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class RvalAlignerCluster {
	public static NeedlemanWunchSolver singleRvalAlignerRun(String seq , RvalSequence Rseq , boolean toPrintAlignment, logObj log) {
		AminoAcidSequence AAseq = new AminoAcidSequence(seq);
		NeedlemanWunchSolver NW = new NeedlemanWunchSolver(Rseq, AAseq, new RvalScoringScheme());
		if (toPrintAlignment) {
			NW.printAlignment();
 			System.out.println();
		}
		log.logString += NW.logAlignment();
		return NW;
	}




	public static void runThread(String swissProtPath, String seqListPath, String profileFilePath) {

		logObj log = new logObj();

//		String swissProtDefaultPath = "/home/zivben/IdeaProjects/Results/uniprot_sprot_10_2015.fasta";
//		String mySeqListDefaultPath = "/home/zivben/IdeaProjects/Results/3JAM_vec3/3JAM_d.fasta";
//		String profileDefaultPath = "/home/zivben/IdeaProjects/Results/3JAM_vec3/3JAM_d_profile.txt";

		String row = "Prot+Chain, Input seq length, *True* length, *True* score, matches, Entries with higher score, Best " +
				"result, \n";



		MySequenceList mySeqList = new MySequenceList(seqListPath); // fasta of single chain

		RvalSequence sourceSeq = new RvalSequence(profileFilePath);	// txt version of CSV "allWithSeq"
		log.logString += "working on profile file: "+profileFilePath+"\n";
		log.logString += mySeqList.fileName().substring(mySeqList.fileName().lastIndexOf("/")+1,mySeqList.fileName()
				.lastIndexOf(".fasta"))+":\t";
		log.logString +="Structure positions:\t" + sourceSeq.size() + "\tLength of *true* seq:\t" + mySeqList.get(0).seq().length() +
				"\n---------------------------------------------------------------------------------------\n";
		System.out.println("\n"+mySeqList.fileName().substring(mySeqList.fileName().lastIndexOf("/")+1,mySeqList.fileName().lastIndexOf(
				".fasta"))+":   Structure positions: " + sourceSeq.size() + "   Length of *true* seq: " + mySeqList
				.get(0).seq().length() + "\n---------------------------------------------------------------------------------------");

		NeedlemanWunchSolver NWsingle = singleRvalAlignerRun(mySeqList.get(0).seq() , sourceSeq, true, log);
		double sameSeqScore = NWsingle.alignmentScore();
		row+= seqListPath.substring(0,seqListPath.indexOf(".fasta"))+",";
		row+= sourceSeq.size()+", "+ mySeqList.get(0).seq().length()+", "+NWsingle.alignmentScore()+", " +NWsingle.bestAlign.matches;

		if (swissProtPath!= null){
			MySequenceList swissProt = new MySequenceList(swissProtPath); // all sequences
			//		MySequenceList swissProt = new MySequenceList("/home/zivben/IdeaProjects/Results/xah.fasta"); // all sequences
			
//			System.out.println("Total num of SwissProt seqs: " + swissProt.size());
//			double largest = 0.0;
//			int indOflargest = -1;
//			int validSeqCounter = 0;
//			int worseScore = 0;
//			log.logString +="\nProcessing SwissProt sequences:\n";
//			for (int ccc=0 ; ccc<swissProt.size() ; ccc++) {
//				if (ccc % 65536 ==0) {System.out.println();}
//				if (ccc % 8192 == 0) {System.out.print(ccc+".. ");}
//				MySequence mySeq = swissProt.get(ccc);
//				if (!(mySeq.seq().indexOf('X')>-1) &&
//						!(mySeq.seq().indexOf('B')>-1) &&
//						!(mySeq.seq().indexOf('J')>-1) &&
//						!(mySeq.seq().indexOf('Z')>-1) &&
//						!(mySeq.seq().indexOf('O')>-1) &&
//						!(mySeq.seq().indexOf('U')>-1)) {
//					if (mySeq.seq().length() >= sourceSeq.size()) {
//						double score = singleRvalAlignerRun(mySeq.seq() , sourceSeq , false, log);
//						if (score>-10000) {
//							validSeqCounter++;
//							if ((score>largest) & (score != sameSeqScore)) {
//								largest = score;
//								indOflargest = ccc;
//							}
//							if (score>sameSeqScore) {
//								worseScore++;
//							}
//						}
//					}
//				}
//			}
//			row += worseScore+", "+largest+", ";
//
//			log.logString +=worseScore + " entries out of valid " + validSeqCounter + " are with better score than the *true*. The " +
//					"best score that is not the *true* sequence is: " + largest+"\n";
//			//		System.out.println(worseScore + " entries out of valid " + validSeqCounter + " are with better score than the *true*. The best score that is not the *true* sequence is: " + largest);
//			if (worseScore>0) {
//				MySequence mySeq = swissProt.get(indOflargest);
//				log.logString +="The best alignment is " + mySeq.title().substring(0,Math.min(mySeq.title().length()-1, 30))+"\n";
//				singleRvalAlignerRun(mySeq.seq() , sourceSeq, true,log);
//				System.out.println("The best alignment is " + mySeq.title().substring(0,Math.min(mySeq.title().length()-1, 30)));
//				//			System.out.println();
//				System.out.println(log.logString);
//			}
//
			
			log.logString += "\nProcessing SwissProt sequences:\n";
//			System.out.println("\nProcessing SwissProt sequences:\n");
			ExecutorCompletionService<Double[]> executorThreads = new ExecutorCompletionService<>(
					new ThreadPoolExecutor(1,3,0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>())
			);
			List<Future<Double[]>> futures = new LinkedList<>();
			int validSwissProt = 0;
			
			for (int ccc = 0; ccc < swissProt.size(); ccc++) {
				//								if (ccc % 65536 == 0) {
				//									setProgress(Math.round((float)ccc / swissProt.size() * 100f));
				//									System.out.print(ccc + ".. ");
				//									publish();
				//								}
				//				if (ccc % 8192 == 0) {
				//					setProgress(Math.round((float)ccc / swissProt.size() * 100f));
				////					System.out.print(ccc + ".. ");
				//					publish();
				//				}
				MySequence mySeq = swissProt.get(ccc);
				if (!(mySeq.seq().indexOf('X') > -1) &&
						!(mySeq.seq().indexOf('B') > -1) &&
						!(mySeq.seq().indexOf('J') > -1) &&
						!(mySeq.seq().indexOf('Z') > -1) &&
						!(mySeq.seq().indexOf('O') > -1) &&
						!(mySeq.seq().indexOf('U') > -1)) {
					if (mySeq.seq().length() >= sourceSeq.size()) {
						final int finalCcc = ccc;
						validSwissProt++;
//						if (validSwissProt % 512 == 0) {
//							System.out.print("\r"+validSwissProt + " valid sequences found");
//						}
						futures.add(executorThreads.submit(() -> {
							AminoAcidSequence targetSeq = new AminoAcidSequence(mySeq.seq());
							NeedlemanWunchSolver NW = new NeedlemanWunchSolver(sourceSeq, targetSeq, new RvalScoringScheme());
							
							return new Double[]{NW.alignmentScore(), (double) finalCcc, NW.getQualityIndex()};
						}));
					}
				}
				
			}
//			System.out.println(validSwissProt + " valid sequences found");
			log.logString +=(validSwissProt + " valid sequences found.\n");
			
			try{
				int counter = 0;
				Future<Double[]> completedThread;
				List<Double[]> resultsList = new LinkedList<>();
				while (counter < validSwissProt) {
					completedThread = executorThreads.take();
					counter++;
					
					resultsList.add(completedThread.get());
//					if (counter % 32 == 0) {
//						System.out.print("\r"+Math.round((float)(counter) / validSwissProt * 100f)+"%\tof sequences checked.");
//					}
				}
				
				
				Collections.sort(resultsList,new ResultsComparator());
				List<Double[]> validList = resultsList.stream().filter(p->p[0]>-1000).collect(Collectors.toList());
				List<Double[]> betterList = validList.stream().filter(FilterPredicates.largerThan(sameSeqScore)).collect(Collectors.toList());
				double largest = validList.get(0)[0];
				int validSeqCounter = validList.size();
				int numBetterScoring = betterList.size();
				
				
				row += numBetterScoring + ", " + largest + ", ";
//				System.out.println();
				log.logString += numBetterScoring + "\tentries out of valid\t" + validSeqCounter + "\tare with better score than the *true*. The " +
						"best score that is not the *true* sequence is:\t" + largest + "\n";
//				System.out.println(numBetterScoring + " entries out of valid " + validSeqCounter + " are with better score than the *true*. The " +
//						"best score that is not the *true* sequence is: " + largest);
//
//				System.out.println("The top 10 best scores are:");
				log.logString += ("The top 10 best scores are:");
				for (int i=0; i<10; i++){
					MySequence mySeq = swissProt.get(validList.get(i)[1].intValue());
					
//					System.out.println("Number " + (i+1) +" is: " + mySeq.title().substring(0, Math.min(mySeq.title().length() - 1, 60))
//							+ "\nWith Quality Index of: " + validList.get(i)[2]);
					
					log.logString += ("Number\t" + (i+1) +"\tis:\t" + mySeq.title().substring(0, Math.min(mySeq.title().length() - 1, 60))
							+ "\nWith Quality Index of:\t" + validList.get(i)[2]+"\n");
					
					singleRvalAlignerRun(mySeq.seq(), sourceSeq, false, log);
//					System.out.println();
					
				}
			
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace();
			}
			
		}
		
		row+="\n";
		File resultCSV = new File(profileFilePath.replace(".txt","_ThreadResults.csv"));
		File logFile = new File(profileFilePath.replace(".txt","_ThreadLog.txt"));

		FileWriter FW;
		try {
			FW = new FileWriter(resultCSV);
			FW.write(row);
			FW.close();

			FW = new FileWriter(logFile);
			FW.write(log.logString );
			FW.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		



	}

	private static class logObj{
		public String logString = "";
		public logObj() {
		}

	}
	
}

class FilterPredicates {
	public static Predicate<Double[]> largerThan(Double trueRes){
		return p -> p[0] > trueRes;
	}
}