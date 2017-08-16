import utils.ExtractMaxValue;
import utils.ScoreUtilities.ScoringGeneralHelpers;
import utils.Scoring.MRC_Score;
import utils.UtilExceptions.MissingChainID;
import utils.molecularElements.SimpleProtein;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


public class Main {

	/**
	 * takes one argument - file with list of chains to process
	 * each line must contain one protein to process, with optional chain designation, in the following format:
	 * <absolute path to PDB file> <abs path to map file> <optional: chain letter, case sensitive>
	 * comment lines in the argument file start with "###" (three hashes).
	 *
	 * @param args path to MRC file
	 */
	public static void main(String[] args) throws MissingChainID {
		// How many threads to run in parrallel (number of available threads multiplied by fraction wanted)
		int cores = (int) Math.ceil(Runtime.getRuntime().availableProcessors() *0.25);
//		int cores = 1;
		System.out.println("Usage: java -jar <jarfile> <argumentFile> <cores> <SCWRL> <swissProtPath>");
		if (args.length > 1) {

			try {
				cores = Integer.parseInt(args[1]);
				ScoringGeneralHelpers.SCWRL_PATH = args[2];
				WorkerThread.swissProtPath = args[3];
				
			} catch (NumberFormatException e) {
				System.err.println("Number of cores: " + args[1] + " must be an integer.");
				System.exit(1);
			}
		}
		Path inputChainList = (new File(args[0])).toPath();
		// create new executor service with a threadpool of the required size.
//		ExecutorService executor = Executors.newFixedThreadPool(cores);
		ExecutorService executor = Executors.newFixedThreadPool(4);


		try {
			// seperate argument file into a list of valid input chains and create new worker thread per chain.
			List<String> chainsToProcess = Files.readAllLines(inputChainList, Charset.defaultCharset());
			for (String chain : chainsToProcess) {
				if (!chain.startsWith("###") && chain.length() >3) {
					String[] tempArgs = chain.split(" ");
					Runnable worker = new WorkerThread(tempArgs);
					executor.execute(worker);
				}
			}
			executor.shutdown();
			while (!executor.isTerminated()) {
			}
			System.out.println("Finished all threads");
		} catch (IOException e) {
			e.printStackTrace();
		}


	}


}

