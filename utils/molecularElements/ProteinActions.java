package utils.molecularElements;

import utils.ScoreUtilities.ScoringGeneralHelpers;
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

import static utils.ScoreUtilities.ScoringGeneralHelpers.*;

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
		SCWRLrunner scwrlRun = new SCWRLrunner(ScoringGeneralHelpers.SCWRL_PATH);

		for (SimpleProtein.ProtChain chain : sourceProtein) {
			for (AminoAcid aminoAcid : chain) {

				for (String newAcid : aAcids) {
					aminoAcid.substituteWith(newAcid);

					// write new processed file
					File fileWithNewRes = new File(outputFolder.getAbsolutePath() + File.separator + sourceProtein
							.getFileName() + "_res_" + aminoAcid.getSeqNum() + "_to_" + newAcid +
							PDB_EXTENSION);

					if (debug) {
						System.out.println("Generating permutation: " + fileWithNewRes.getName());
					}
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

		if (name.equals("ALA") || name.trim().equals("A"))
			return 0;
		if (name.equals("ARG") || name.trim().equals("R"))
			return 1;
		if (name.equals("ASN") || name.trim().equals("N"))
			return 2;
		if (name.equals("ASP") || name.trim().equals("D"))
			return 3;
		if (name.equals("CYS") || name.trim().equals("C"))
			return 4;
		if (name.equals("GLU") || name.trim().equals("E"))
			return 5;
		if (name.equals("GLN") || name.trim().equals("Q"))
			return 6;
		if (name.equals("GLY") || name.trim().equals("G"))
			return 7;
		if (name.equals("HIS") || name.trim().equals("H"))
			return 8;
		if (name.equals("ILE") || name.trim().equals("I"))
			return 9;
		if (name.equals("LEU") || name.trim().equals("L"))
			return 10;
		if (name.equals("LYS") || name.trim().equals("K"))
			return 11;
		if (name.equals("MET") || name.trim().equals("M"))
			return 12;
		if (name.equals("PHE") || name.trim().equals("F"))
			return 13;
		if (name.equals("PRO") || name.trim().equals("P"))
			return 14;
		if (name.equals("SER") || name.trim().equals("S"))
			return 15;
		if (name.equals("THR") || name.trim().equals("T"))
			return 16;
		if (name.equals("TRP") || name.trim().equals("W"))
			return 17;
		if (name.equals("TYR") || name.trim().equals("Y"))
			return 18;
		if (name.equals("VAL") || name.trim().equals("V"))
			return 19;

		if (name.equals("SEC") || name.trim().equals("U"))
			return 4;
			//TODO - fix this 21st amino acid thing.

		else
			throw new InvalidPropertiesFormatException("Bad AminoAcid Name: " + name);
	}
}