package utils.scwrlIntegration;

import java.io.File;
import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import static utils.fileUtilities.FileProcessor.*;

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

		List<File> fileNames = new ArrayList<>();
		try (DirectoryStream<Path> directoryStream = Files.newDirectoryStream(tempFolder.toPath())) {
			for (Path path : directoryStream) {
				fileNames.add(path.toFile());
			}
		} catch (IOException ex) {
			throw ex;
		}

		File targetFolder = makeFolder(tempFolder, "_SCWRL");
		for (File fileName : fileNames) {
			scwrlRunOnce(fileName, targetFolder);
		}

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
