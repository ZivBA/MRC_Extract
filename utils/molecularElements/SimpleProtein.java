package utils.molecularElements;

import utils.fileUtilities.FileProcessor;
import utils.scwrlIntegration.SCWRLactions;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static utils.fileUtilities.FileProcessor.*;

/**
 * Created by Ziv_BA on 30/07/2015.
 */
public class SimpleProtein implements Iterable<SimpleProtein.ProtChain> {

	public File processingFolder;
	private File source;
	private String fileName;
	private int sequenceBias;
	private List<ProtChain> protChains; // list of the seperate AminoAcid chains
	private List<String> hetAtmAndFooter;   // array of the remaining HeteroAtoms + footer tags.
	private double[][] intensityMatrix;

	/**
	 * constructor for SimpleProtein from PDB file.</br>
	 * creates atoms for each ATOM line, combine atoms by resSeq as residues, all residues of a specific
	 * chain in a ProtChain object. each chain ends with a TER line. </br>
	 * hetero atoms, master tag and end tag are all injected into the hetAtmAndFooter String array without
	 * processing.
	 *
	 * @param pdbFile File Object pointing to a PDB file.
	 * @throws IOException if an error occurs while reading the file.
	 */
	public SimpleProtein(File pdbFile) throws IOException {
		source = pdbFile;
		fileName = pdbFile.getName().substring(0, pdbFile.getName().indexOf(PDB_EXTENSION));

		protChains = new ArrayList<>();
		// create a chain 2D array and the first chain list.
		List<ArrayList<String>> chains = new ArrayList<>();
		chains.add(new ArrayList<String>());

		List<String> hetAtm = new ArrayList<>();
		int chainCounter = 0;

		List<String> PDBStrArr = Files.readAllLines(pdbFile.toPath(), Charset.defaultCharset());


		// iterate over the entire PDB file, populate the chain lists with ATOM lines, and hetAtom with
		// HETATM lines.
		// for each TER entry, add a new list for that chain.
		for (int i = 0; i < PDBStrArr.size(); i++) {
			String lineInFile = PDBStrArr.get(i);
			if (lineInFile.startsWith("ATOM")) {
				chains.get(chainCounter).add(lineInFile);
			} else if (lineInFile.startsWith("TER")) {
				chains.get(chainCounter).add(lineInFile);
				chainCounter++;
				chains.add(new ArrayList<String>());
			} else if (lineInFile.matches("(" + FOOTER_TAGS + ").*")) {
				hetAtm.add(lineInFile);
			}

		}

		for (List<String> chain : chains) {
			if (chain != null && chain.size() > 1) {
				protChains.add(new ProtChain(chain));
			}
		}

		hetAtmAndFooter = hetAtm;
		sequenceBias = protChains.get(0).residues.get(0).getSeqNum();


	}

	public void setIntensityMatrix(double[][] intensityMatrix) {
		this.intensityMatrix = intensityMatrix;
	}

	public String getFileName() {
		return fileName;
	}

	public int getSequenceBias() {
		return sequenceBias;
	}

	public File getSource() {
		return source;
	}

	/**
	 * dump all the atom sequences and footer back to a PDB file.
	 *
	 * @param destination
	 * @throws IOException
	 */
	public void writePDB(File destination) throws IOException {
		FileWriter FW = new FileWriter(destination);

		// write structural atoms
		for (ProtChain chain : protChains) {
			for (AminoAcid acid : chain) {
				for (SimpleAtom atom : acid) {
					FW.write(atom.getOriginalString() + "\n");
				}
			}

		}
		//write HeteroAtoms and footer tags.
		for (String line : hetAtmAndFooter) {
			FW.write(line);
		}
		FW.close();
	}

	public void createPermutations() throws IOException {
		processingFolder = FileProcessor.makeFolder(source, "_temp");
		ProteinActions.iterateAllAcidsToFile(this, processingFolder);
		SCWRLactions.genSCWRLforFolder(processingFolder);

	}

	public void createPermutations(File tempFolder) throws IOException {
		this.processingFolder = tempFolder;
		makeFolder(processingFolder);

		ProteinActions.iterateAllAcidsToFile(this, processingFolder);
		SCWRLactions.genSCWRLforFolder(processingFolder);


	}

	public File getProcessingFolder() {
		return processingFolder;
	}

	public void setProcessingFolder(File processingFolder) {
		this.processingFolder = processingFolder;
	}

	@Override
	public Iterator<ProtChain> iterator() {
		return protChains.iterator();
	}

	public int getLegnth() {
		int protLength = 0;
		for (ProtChain chain : protChains) {
			protLength += chain.residues.size();
		}
		return protLength;

	}

	/**
	 * helper class to bulk residues together in respective chains.
	 * also performs the actual processing from string array to molecular elements.
	 */
	public class ProtChain implements Iterable<AminoAcid> {
		private char chainID;
		private List<AminoAcid> residues = new ArrayList<>();

		/**
		 * constructor creating a chain from a list of strings (assume all strings are for a single chain)
		 *
		 * @param sourceList
		 */
		protected ProtChain(List<String> sourceList) {
			List<String> tempAtomList = new ArrayList<>();
			int workingResSeq = Integer.valueOf(sourceList.get(0).substring(RES_SEQ_START, RES_SEQ_END + 1)
					.trim());

			for (int i = 0; i < sourceList.size(); i++) {

				String tempCurrLine = sourceList.get(i);
				int tempCurrResSeq = Integer.valueOf(
						tempCurrLine.substring(RES_SEQ_START, RES_SEQ_END + 1).trim());

				if (tempCurrResSeq == workingResSeq) {
					tempAtomList.add(tempCurrLine);
				} else {
					workingResSeq = tempCurrResSeq;
					residues.add(new AminoAcid(tempAtomList));
					tempAtomList.clear();
					tempAtomList.add(tempCurrLine);
				}

			}

			chainID = residues.get(0).getChainID();

		}

		public char getChainID() {
			return chainID;
		}

		protected void addRes(AminoAcid res) {
			residues.add(res);
		}

		@Override
		public Iterator<AminoAcid> iterator() {
			return residues.iterator();
		}
	}
}
