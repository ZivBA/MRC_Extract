package utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import static utils.FileProcessor.*;

/**
 * Created by zivben on 04/08/15.
 */
public class ProteinActions {

	/**
	 * strips a PDB file from all residue atoms except basic backbone (Ca, C, N, O) and replaces all AA to ALA
	 * @param source
	 * @param dest
	 * @throws IOException
	 */
	public static void stripAndAllALA(File source, File dest) throws IOException {


		List<String> PDBStrArr = Files.readAllLines(source.toPath(), Charset.defaultCharset());
		Iterator<String> PDBIter = PDBStrArr.iterator();
		Pattern acceptableNames = Pattern.compile("\\s*("+ALLOWED_ATOMS+")\\s*");
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

	public static void iterateAllAcids(SimpleProtein prot, File tempFolder) throws IOException {
		String srcName = prot.getSource().getName();



		// go over all AA in the protein
		for (SimpleProtein.ProtChain chain : prot){
			for (AminoAcid aminoAcid : chain) {

				// for each Amino acid, replace with every possible AA
				for (String acid : aAcids){
					aminoAcid.substituteWith(acid);
					prot.writePDB(new File(tempFolder.getAbsolutePath()+File.separator + srcName +
							"_res_"+aminoAcid.getSeqNum()+"_to_"+acid));

					//reset to ALA
					aminoAcid.substituteWith("ALA");

				}
			}

		}
	}

	public static void genSCWRLforFolder(File tempFolder) {

	}
}