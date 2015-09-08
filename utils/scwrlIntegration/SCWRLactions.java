package utils.scwrlIntegration;

import java.io.File;
import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import static utils.ScoreUtilities.ScoringGeneralHelpers.*;

/**
 * Created by zivben on 06/08/15.
 */
public class SCWRLactions {

	/**
	 * iterate all files in folder and scwrl them.
	 *
	 * @param tempFolder
	 */
	public static void genSCWRLforFolder(File tempFolder) throws IOException {

		System.out.println("******************************************************************");
		System.out.println("Processing SCWRL input folder: \n" + tempFolder.getAbsolutePath());
		System.out.println("******************************************************************");
		List<File> fileNames = new ArrayList<>();
		long startTime = System.currentTimeMillis();
		float tempTime;
		try (DirectoryStream<Path> directoryStream = Files.newDirectoryStream(tempFolder.toPath())) {
			for (Path path : directoryStream) {
				if (!path.getFileName().toString().endsWith("_SCWRLed.pdb")) {
					fileNames.add(path.toFile());
				}
			}
		} catch (IOException ex) {
			throw ex;
		}

		int numOfFiles = fileNames.size();
		int blockSize = (int) (numOfFiles / (Math.ceil(Math.log10(numOfFiles))));
		int filesScwrled = 0;
		SCWRLrunner scwrl = new SCWRLrunner(SCWRL_PATH);
		for (File fileName : fileNames) {
			File SCWRLFile = new File(fileName.getAbsolutePath().replace(".pdb", "_SCWRLed.pdb"));
			if (!SCWRLFile.exists()) {
				scwrl.runScwrl(fileName, SCWRLFile);
			}
			fileName.delete();
			filesScwrled++;
			if (filesScwrled % blockSize == 0) {
				tempTime = System.currentTimeMillis();
				float elapsed = (tempTime - startTime) / 1000f;
				System.out.println("Processed " + filesScwrled + " files out of " + numOfFiles + "\nthis batch took: "
						+ String.valueOf(elapsed));
				System.out.println("Should probably take about " + ((numOfFiles - filesScwrled) / blockSize *
						elapsed) +
						" Seconds");
			}
		}
		long stopTime = System.currentTimeMillis();
		float elapsedTime = (stopTime - startTime) / 1000f;
		System.out.println("******************************************************************");
		System.out.println("Generated: " + fileNames.size() + " Files in: " + elapsedTime + " seconds");
		System.out.println("******************************************************************");
		System.out.println("SCWRL execution terminated!");


	}

	public static void scwrlRunOnce(File inputFile, File targetFolder) throws IOException {

		SCWRLrunner scwrl = new SCWRLrunner(SCWRL_PATH);
		File newScwrlFile;


		if (debug) {

			newScwrlFile = new File(
					targetFolder.getAbsolutePath() + File.separator + inputFile.getName().replaceFirst(
							PDB_EXTENSION + "+$", "_SCWRL" + PDB_EXTENSION));

		} else {
			newScwrlFile = inputFile;
		}


		scwrl.runScwrl(inputFile, newScwrlFile);

	}


}
