package utils.fileUtilities;


import utils.ExtractMaxValue;

import java.io.File;
import java.io.IOException;

/**
 * Created by Ziv_BA on 30/07/2015.
 */
public class FileProcessor {
	//path to SCWRL
	public static final String SCWRL_PATH = "/home/zivben/SCWRL4/Scwrl4";
	/**
	 * constatns from PDB format. (positions are "-1" since format starts at '1' and computers start at '0'
	 */
	public static final int ATOM_NUM_START = 6, ATOM_NUM_END = 10, ATOM_NAME_START = 12, ATOM_NAME_END = 15,
			RES_NAME_START = 17, RES_NAME_END = 19, CHAIN_ID = 21, RES_SEQ_START = 22, RES_SEQ_END = 25;
	public static final String FOOTER_TAGS = "HETATM|MASTER|END", ALLOWED_ATOMS = "N|CA|C|O", PDB_EXTENSION
			= ".pdb";
	//consts
	// keep ALA as first item, as this is being called explicitly by stripAllRes method.
	public static final String[] aAcids = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
			"ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
	public static boolean debug = false;
	private File source;
	private File dest;

	/**
	 * Constructor - gets File obj, sets destination.
	 *
	 * @param toProcess
	 */
	public FileProcessor(File toProcess, boolean debug) {
		FileProcessor.debug = debug;
		source = toProcess;
		dest = new File(source.getAbsolutePath().replaceFirst("[.][pdb]+$", "_stripped" + PDB_EXTENSION));
	}

	public static File makeFolder(File source, String newFolder) throws IOException {

		File requestedFolder = new File(source.getParent() + File.separator + newFolder);

		if (requestedFolder.isDirectory()) {
			System.out.println("Requested folder already exists at:\n" + requestedFolder.getAbsolutePath());
		} else {
			if (requestedFolder.mkdir()) {
				System.out.println("Created requested folder \n" + requestedFolder.getAbsolutePath());
			} else {
				System.out.println("Requested Folder not created");
				throw new IOException();
			}
		}

		return requestedFolder;
	}

	public static File makeFolder(File requestedFolder) throws IOException {

		if (requestedFolder.isDirectory()) {
			System.out.println("Requested folder already exists at:\n" + requestedFolder.getAbsolutePath());
		} else {
			if (requestedFolder.mkdir()) {
				System.out.println("Created requested folder \n" + requestedFolder.getAbsolutePath());
			} else {
				System.out.println("Requested Folder not created");
				throw new IOException();
			}
		}

		return requestedFolder;
	}

	/**
	 * strips all residue atoms from CBeta onwards, replace all residue names with ALA;
	 *
	 * @throws IOException
	 */


	public float[][] scoreProteinAgainstMRC(MRC_Map_New map) {


		return null;
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
		System.out.println("heighest intensity was: " + results[0] + " was detected at: " + results[1] +
				"," +
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
}