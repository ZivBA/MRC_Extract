package utils.molecularElements;

import utils.fileUtilities.FileProcessor;
import utils.scwrlIntegration.SCWRLrunner;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.InvalidPropertiesFormatException;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static utils.fileUtilities.FileProcessor.*;

/**
 * Created by zivben on 04/08/15.
 */
public class ProteinActions {

	public static void stripAndAllALAToObject(SimpleProtein sourceProtein) {

		for (SimpleProtein.ProtChain chain : sourceProtein) {
			for (AminoAcid acid : chain) {
				acid.substituteWith("ALA");
				acid.strip();
			}

		}

	}

	/**
	 * strips a PDB file from all residue atoms except basic backbone (Ca, C, N, O) and replaces all AA to
	 * ALA. saves to file for all to enjoy.
	 *
	 * @param source
	 * @param dest
	 * @throws IOException
	 */
	public static void stripAndAllALAToFile(File source, File dest) throws IOException {


		List<String> PDBStrArr = Files.readAllLines(source.toPath(), Charset.defaultCharset());
		Iterator<String> PDBIter = PDBStrArr.iterator();
		Pattern acceptableNames = Pattern.compile("\\s*(" + ALLOWED_ATOMS + ")\\s*");
		Matcher matchNames = acceptableNames.matcher("");
		FileWriter writer = new FileWriter(dest);
		while (PDBIter.hasNext()) {
			String line = PDBIter.next() + "\n";
			if (line.startsWith("ATOM")) {
				if (matchNames.reset(line.substring(ATOM_NAME_START, ATOM_NAME_END + 1)).matches()) {
					line = line.substring(0, RES_NAME_START) + aAcids[0] + line.substring(RES_NAME_END + 1);
					writer.write(line);
				}
			} else {
				writer.write(line);
			}


		}
		writer.close();
	}

	public static File iterateAndScwrl(SimpleProtein sourceProtein) throws IOException {


		File outputFolder = makeSubFolderAt(sourceProtein.getSource(), "_scwrlFiles");
		SCWRLrunner scwrlRun = new SCWRLrunner(FileProcessor.SCWRL_PATH);

		for (SimpleProtein.ProtChain chain : sourceProtein) {
			for (AminoAcid aminoAcid : chain) {

				for (String newAcid : aAcids) {
					aminoAcid.substituteWith(newAcid);
					// write new processed file
					File fileWithNewRes = new File(outputFolder.getAbsolutePath() + File.separator + sourceProtein
							.getFileName() + "_res_" + aminoAcid.getSeqNum() + "_to_" + newAcid +
							PDB_EXTENSION);

					sourceProtein.writePDB(fileWithNewRes);
					scwrlRun.runScwrl(fileWithNewRes, fileWithNewRes);


					//reset to ALA
					aminoAcid.substituteWith("ALA");

				}


			}
		}


		return outputFolder;

	}


	/**
	 * take a protein and start iterating all over it's residues. at each position create a new file with
	 * that residue replaced to each of the 20 AA.
	 *
	 * @param prot         input protein
	 * @param outputFolder output folder
	 * @throws IOException
	 */
	public static void iterateAllAcidsToFile(SimpleProtein prot, File outputFolder) throws IOException {

		// go over all AA in the protein
		for (SimpleProtein.ProtChain chain : prot) {
			for (AminoAcid aminoAcid : chain) {

				// for each Amino acid, replace with every possible AA
				for (String acid : aAcids) {
					aminoAcid.substituteWith(acid);
					// write new processed file
					File fileWithNewRes = new File(outputFolder.getAbsolutePath() + File.separator + prot
							.getFileName() + "_res_" + aminoAcid.getSeqNum() + "_to_" + acid +
							PDB_EXTENSION);

					prot.writePDB(fileWithNewRes);

					//reset to ALA
					aminoAcid.substituteWith("ALA");

				}
			}

		}
	}

	public static int acidToIndex(String name) throws InvalidPropertiesFormatException {

		if (name.equals("ALA"))
			return 0;
		if (name.equals("ARG"))
			return 1;
		if (name.equals("ASN"))
			return 2;
		if (name.equals("ASP"))
			return 3;
		if (name.equals("CYS"))
			return 4;
		if (name.equals("GLU"))
			return 5;
		if (name.equals("GLN"))
			return 6;
		if (name.equals("GLY"))
			return 7;
		if (name.equals("HIS"))
			return 8;
		if (name.equals("ILE"))
			return 9;
		if (name.equals("LEU"))
			return 10;
		if (name.equals("LYS"))
			return 11;
		if (name.equals("MET"))
			return 12;
		if (name.equals("PHE"))
			return 13;
		if (name.equals("PRO"))
			return 14;
		if (name.equals("SER"))
			return 15;
		if (name.equals("THR"))
			return 16;
		if (name.equals("TRP"))
			return 17;
		if (name.equals("TYR"))
			return 18;
		if (name.equals("VAL"))
			return 19;

		else
			throw new InvalidPropertiesFormatException("Bad AminoAcid Name");
	}
}