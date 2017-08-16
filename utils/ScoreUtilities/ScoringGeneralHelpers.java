package utils.ScoreUtilities;


import utils.ExtractMaxValue;
import utils.molecularElements.SimpleProtein;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by Ziv_BA on 30/07/2015.
 */
public class ScoringGeneralHelpers {
	//path to SCWRL
	public static String SCWRL_PATH = "/home/zivben/SCWRL4/Scwrl4";
	/**
	 * constatns from PDB format. (positions are "-1" since format starts at '1' and computers start at '0'
	 */
	public static final int ATOM_NUM_START = 6, ATOM_NUM_END = 10, ATOM_NAME_START = 12, ATOM_NAME_END = 15,
			RES_NAME_START = 17, RES_NAME_END = 19, CHAIN_ID = 21, RES_SEQ_START = 22, RES_SEQ_END = 25;
	public static final String FOOTER_TAGS = "HETATM|MASTER|END", ALLOWED_ATOMS = "N|CA|C|O", PDB_EXTENSION
			= ".pdb";
	//consts
	// keep ALA as first item, as this is being called explicitly by stripAllRes method.
	public static final String[] aAcids = {"ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE",
			"LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"};

	public static final char[] singleLetters = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
			'R', 'S', 'T', 'V', 'W', 'Y'};

	// vector designating which amino acids have negative signal. those are to be multiplied by "-1" to normalize the
	// zvalues for further processing. vector is 20 positions in order of AA same as seen above.
//	public static final double[] vector2 = {1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0,
//			1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0,};
//
//	public static final double[] vector3 = {1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0,
//			1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0,};
//
//	public static final double[] vector4 = {-1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
//			1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0,};


	public static double[] latestVector = {-1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0,};


	public static double[] normalVector = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,};


	public static boolean debug = false;
	private File source;
	private File dest;

	/**
	 * Constructor - gets File obj, sets destination.
	 *
	 * @param toProcess
	 */
	public ScoringGeneralHelpers(File toProcess, boolean debug) {
		ScoringGeneralHelpers.debug = debug;
		source = toProcess;
		dest = new File(source.getAbsolutePath().replaceFirst("[.][pdb]+$", "_stripped" + PDB_EXTENSION));
	}

	public static File makeSubFolderAt(File sourceFile, String targetSubFolder) throws IOException {
		File requestedFolder;
		if (sourceFile.isFile()) {
			requestedFolder = new File(sourceFile.getParent() + File.separator + targetSubFolder);
		} else if (sourceFile.isDirectory()) {
			requestedFolder = new File(sourceFile.getAbsolutePath() + File.separator + targetSubFolder);
		} else {
			throw new FileNotFoundException("source folder is not a proper path.");
		}

		if (requestedFolder.isDirectory()) {
			if (debug)
				System.out.println("Requested folder already exists at:\n" + requestedFolder.getAbsolutePath());
		} else {
			if (requestedFolder.mkdir()) {
				if (debug)
					System.out.println("Created requested folder \n" + requestedFolder.getAbsolutePath());
			} else {
				if (debug)
					System.out.println("Requested Folder not created");
				throw new IOException();
			}
		}

		return requestedFolder;
	}

	public static File makeFolder(File requestedFolder) throws IOException {

		if (requestedFolder.isDirectory()) {
			if (debug)
				System.out.println("Requested folder already exists at:\n" + requestedFolder.getAbsolutePath());
		} else {
			if (requestedFolder.mkdir()) {
				if (debug)
					System.out.println("Created requested folder \n" + requestedFolder.getAbsolutePath());
			} else {
				if (debug)
					System.out.println("Requested Folder not created");
				throw new IOException();
			}
		}

		return requestedFolder;
	}

	public static void multiplyMatrixByVector(double[][] allZvalueMatrix, double[] vector) {
		for (int i = 0; i < allZvalueMatrix.length; i++) {
			for (int j = 0; j < allZvalueMatrix[i].length; j++) {
				allZvalueMatrix[i][j] *= vector[i];
			}
		}
	}

	public static void multiplyMatrixByVector(SimpleProtein.ProtChain chain, double[] vector) {
		for (int i = 0; i < chain.allZvalueMatrix.length; i++) {
			for (int j = 0; j < chain.allZvalueMatrix[i].length; j++) {
				chain.allZvalueMatrix[i][j] *= vector[i];
			}
		}
		for (int i = 0; i < chain.trueZvalues.length; i++) {
			chain.trueZvalues[i] *= vector[chain.originalPositions[i]];
		}
	}

	public static void multiplyMatrixByVector(double[] truzvaltmp, double[] vector) {
		for (int i = 0; i < truzvaltmp.length; i++) {
			truzvaltmp[i] *= vector[i];
		}
	}

	private void getMaxIntensityValueFromMRCFile(String arg) {
		// create MRC object with helper class MRC_map
		MRC_Map_New mrcMap = null;
		try {
			mrcMap = new MRC_Map_New(arg);
		} catch (RuntimeException e) {
			System.out.println("Problem with input file, error trace follows:");
			System.out.println(e.getMessage() + "\n");
			e.printStackTrace();
		}

		float[] results = ExtractMaxValue.getMaxValue(mrcMap);
		System.out.println("heighest intensity was: " + results[0] + " was detected at: " + results[1] + "," +
				results[2] + "," + results[3]);

		if (mrcMap != null) {
			System.out.println(mrcMap.val(results[1], results[2], results[3]));
		} else {
			System.out.println("Something wrong with MRC_map object.");
		}
	}

	public File getDest() {
		return dest;
	}

	public File getSource() {
		return source;
	}


	public static double[][] csvToMatrix(File resultCSV) {
		double[][] resultMatrix;
		try {
			List<String> csvList = Files.readAllLines(resultCSV.toPath(), Charset.defaultCharset());
			resultMatrix = new double[csvList.get(0).split(",").length][csvList.size()];
			for (int i = 0; i < csvList.size(); i++) {
				String[] tempLine = csvList.get(i).split(",");
				for (int j = 0; j < csvList.get(0).split(",").length; j++) {
					resultMatrix[j][i] = Double.parseDouble(tempLine[j]);
				}
			}

		} catch (IOException e) {
			if (debug)
				System.out.println("Cannot read requested CSV");
			return null;
		}
		return resultMatrix;
	}
}