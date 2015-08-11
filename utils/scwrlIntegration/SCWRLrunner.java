package utils.scwrlIntegration;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by zivben on 06/08/15.
 */
public class SCWRLrunner {

	File scwrlExe;

	public SCWRLrunner(String pathToScwrlExe) throws FileNotFoundException {

		scwrlExe = new File(pathToScwrlExe);

		if (!scwrlExe.isFile()) {
			throw new FileNotFoundException(pathToScwrlExe);
		}

	}

	public String[] runScwrl(File input, File output) throws IOException {

		if (!output.isFile()) {
			output.createNewFile();
			output.setWritable(true);
		}


		Process process = Runtime.getRuntime().exec(scwrlExe.getAbsolutePath() +
				" -i " + input.getAbsolutePath() +
				" -o " + output.getAbsolutePath() +
				" -h");

		BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()));

		List stdOutput = new ArrayList();
		String line;
		while ((line = br.readLine()) != null) {
			stdOutput.add(line);
		}
		return (String[]) stdOutput.toArray(new String[stdOutput.size()]);
	}
}
