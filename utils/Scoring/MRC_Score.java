package utils.Scoring;

import utils.ScoreUtilities.MRC_Map_New;
import utils.ScoreUtilities.ScoringGeneralHelpers;
import utils.UtilExceptions.MissingChainID;
import utils.molecularElements.AminoAcid;
import utils.molecularElements.ProteinActions;
import utils.molecularElements.SimpleAtom;
import utils.molecularElements.SimpleProtein;
import utils.scwrlIntegration.SCWRLactions;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.InvalidPropertiesFormatException;
import java.util.List;

import static org.apache.commons.lang3.StringUtils.substringBetween;
import static utils.molecularElements.ProteinActions.acidToIndex;

/**
 * Created by zivben on 09/08/15.
 */
public class MRC_Score {
	List<Integer[]> originalPos;
	private MRC_Map_New myMap;
	private SimpleProtein myProt;
	private char requestedChain;
	private double[][] proteinIntensityValueMatrix;

	public MRC_Score(MRC_Map_New myMap, SimpleProtein myProt) {
		this.myMap = myMap;
		this.myProt = myProt;
		proteinIntensityValueMatrix = new double[20][myProt.getNumChains()];

	}

	public MRC_Score(String mapPath, String protPath) throws IOException {
		this(new MRC_Map_New(mapPath), new SimpleProtein(new File(protPath)));
	}

	public MRC_Score(String mapPath, String protPath, String requestedChain) throws IOException {
		this(new MRC_Map_New(mapPath), new SimpleProtein(new File(protPath), requestedChain));
		this.requestedChain = requestedChain.charAt(0);

	}

	public MRC_Map_New getMyMap() {
		return myMap;
	}

	//	public static MRC_Score StartFromScratch(ScoringGeneralHelpers FP, String mrcpath) {
	//
	//		try {
	//			ProteinActions.stripAndAllALAToFile(FP.getSource(), FP.getDest());
	//			SimpleProtein processedProt = new SimpleProtein(FP.getDest());
	//
	//			MRC_Score scoreMap = new MRC_Score(new MRC_Map_New(mrcpath), processedProt);
	//			scoreMap.scoreProtein();
	//			scoreMap.calcZvalue();
	//			scoreMap.createCSVs();
	//			return scoreMap;
	//		} catch (IOException e) {
	//			e.printStackTrace();
	//			return null;
	//		}
	//
	//	}

	/**
	 * process the myProt SimpleProtein object, strip residues and iterate all permutations of amino acids.
	 * run SCWRL external utility for each permutation and return the intensity value matrix result.
	 *
	 * @return 2D array of double precision values representing intensity values for each amino acid at each position.
	 * @throws IOException
	 */
	public double[][] scoreProtein() throws IOException, MissingChainID {
		System.out.println("Starting protein scoring, saving original positions.");
		myProt.saveOriginalPositions();
		originalPos = myProt.getOriginalPositions();

		System.out.println("Stripping all amino acid resiudes and setting to ALA.");
		ProteinActions.stripAndAllALAToObject(myProt);
		System.out.println("Iterating all acid permutations and creating SCWRL input files");
		File scwrlOutput = ProteinActions.iterateAcids(myProt);

		File chainSubFolders[] = scwrlOutput.listFiles();
		for (File chain : chainSubFolders) {
			if (chain.isDirectory()) {
				if (chain.getAbsolutePath().endsWith(File.separator + requestedChain) || requestedChain == '\0') {
					SCWRLactions.genSCWRLforFolder(chain);
				}
			}
		}

		processSCWRLfolder(scwrlOutput, requestedChain);

		return proteinIntensityValueMatrix;
	}

	private void processSCWRLfolder(File processingFolder, char requestedChain) throws IOException, MissingChainID {

		SimpleProtein tempProt;

		List<File> chainFolders = new ArrayList<>();
		try (DirectoryStream<Path> directoryStream = Files.newDirectoryStream(processingFolder.toPath())) {
			for (Path path : directoryStream) {
				if (path.toString().endsWith(File.separator + requestedChain)) {
					chainFolders.add(path.toFile());
				} else if (requestedChain == '\0') {
					chainFolders.add(path.toFile());
				}
			}
		} catch (IOException ex) {
			throw ex;
		}

		for (File chain : chainFolders) {
			try (DirectoryStream<Path> directoryStream = Files.newDirectoryStream(chain.toPath())) {
				for (Path path : directoryStream) {
					tempProt = new SimpleProtein(path.toFile());
					scoreSingleScwrl(tempProt);
				}
			} catch (IOException | MissingChainID ex) {
				throw ex;
			}


		}
	}

	private void processSCWRLfolder(File processingFolder) throws IOException, MissingChainID {

		SimpleProtein tempProt;

		List<File> chainFolders = new ArrayList<>();
		try (DirectoryStream<Path> directoryStream = Files.newDirectoryStream(processingFolder.toPath())) {
			for (Path path : directoryStream) {
				chainFolders.add(path.toFile());
			}
		} catch (IOException ex) {
			throw ex;
		}

		for (File chain : chainFolders) {
			try (DirectoryStream<Path> directoryStream = Files.newDirectoryStream(chain.toPath())) {
				for (Path path : directoryStream) {
					tempProt = new SimpleProtein(path.toFile());
					scoreSingleScwrl(tempProt);
				}
			} catch (IOException | MissingChainID ex) {
				throw ex;
			}


		}
	}

	private void scoreSingleScwrl(SimpleProtein tempProt) throws InvalidPropertiesFormatException, MissingChainID {


		int testResiduePosition = Integer.valueOf(substringBetween(tempProt.getFileName(), "_res_", "_"));
		int testResidueIndex = acidToIndex(
				substringBetween(tempProt.getFileName(), "_to_", "_SCWRLed"));

		for (SimpleProtein.ProtChain chain : tempProt) {
			if (chain.getChainID() == requestedChain || requestedChain == '\0') {
				for (AminoAcid residue : chain) {

					if (residue.getSeqNum() == testResiduePosition && acidToIndex(
							residue.getName()) == testResidueIndex) {
						double resSum = 0;
						double backBoneSum = 0;
						for (SimpleAtom atom : residue) {
							if (atom.isBackbone()) {
								backBoneSum += scoreSingleAtom(atom);
							} else {
								resSum += scoreSingleAtom(atom);
							}
						}

						SimpleProtein.ProtChain originalChain = myProt.getChain(residue.getChainID());
						if (originalChain == null)
							throw new MissingChainID(residue.getChainID(), tempProt);
						originalChain.intensityValueMatrix[acidToIndex(
								residue.getName())][residue.getSeqNum() - originalChain.getSequenceBias()] =
							/*backBoneSum  +*/ resSum;


					}
				}
			}
		}


	}
	
	private double scoreSingleAtom(SimpleAtom atom) {
		float[] coords = atom.getAtomCoords();

		try {
			atom.setAtomScore(myMap.val(coords[0], coords[1], coords[2]));
		} catch (RuntimeException e) {
			System.out.println("Warning, PDB contains coordinate value outside MRC map scope.");
			System.out.println("Atom " + atom.getName() + " from residue " + atom.getaAcidName() + " at coordinates " +
					Arrays.toString(atom.getAtomCoords()) + " is the problem.");
			System.out.println("problematic line from PDB is: \n" + atom.getOriginalString() + "\n");
			atom.setAtomScore(0.0);
		}
		return atom.getAtomScore();
	}

	public void calcZvalue() throws InvalidPropertiesFormatException {
		System.out.println("Calculating Z-Values");
		for (SimpleProtein.ProtChain chain : myProt) {
			zValueHelper(chain);
		}

	}

	private void zValueHelper(SimpleProtein.ProtChain chain) {
		double tempAvg[] = new double[20];
		double tempStD[] = new double[20];
		chain.zValueMatrix = new double[chain.intensityValueMatrix.length][chain.intensityValueMatrix[0].length];


		// calc tempAvg per column in the intensity value matrix ( avarage of scores per amino acid in every pos)
		for (int i = 0; i < chain.intensityValueMatrix.length; i++) {
			for (int j = 0; j < chain.intensityValueMatrix[i].length; j++) {
				tempAvg[i] += chain.intensityValueMatrix[i][j];
			}
			tempAvg[i] = tempAvg[i] / chain.intensityValueMatrix[i].length;

			//calc standard deviation for each column
			for (int j = 0; j < chain.intensityValueMatrix[i].length; j++) {
				tempStD[i] += Math.pow(chain.intensityValueMatrix[i][j] - tempAvg[i], 2);
			}
			tempStD[i] = Math.sqrt(tempStD[i] / chain.intensityValueMatrix[i].length);
		}

		// calc Z-Value for each discrete acid in every position
		for (int i = 0; i < chain.intensityValueMatrix.length; i++) {
			for (int j = 0; j < chain.intensityValueMatrix[i].length; j++) {

				if (chain.originalPositions[j] == i) {
					chain.originalAcidsScore[j] = (chain.intensityValueMatrix[i][j] - tempAvg[i]) / tempStD[i];
					chain.zValueMatrix[i][j] = 0;
				} else {
					chain.zValueMatrix[i][j] = (chain.intensityValueMatrix[i][j] - tempAvg[i]) / tempStD[i];
				}

			}
		}

		Integer[] suspectedCorrectPositions = chain.getOriginalPositions();

		for (int i = 0; i < suspectedCorrectPositions.length; i++) {
			chain.zValueCorrect[i] = chain.zValueMatrix[suspectedCorrectPositions[i]][i];
		}
	}

	public void createCSVs() throws IOException {

		System.out.println("Creating CSV Files");
		File tempCSVfolder = ScoringGeneralHelpers.makeFolder(new File(myProt.getSource().getParent() + File
				.separator + "tempCSVs"));


		for (SimpleProtein.ProtChain chain : myProt) {

			File resultCSV = new File(tempCSVfolder.getAbsolutePath() + File.separator + myProt.getFileName()
					+ "_Chain_" + chain.getChainID() + "_resultMatrix.csv");
			File zscoreCSV = new File(tempCSVfolder.getAbsolutePath() + File.separator + myProt.getFileName()
					+ "_Chain_" + chain.getChainID() + "_zscore.csv");
			File zscoreCorrect = new File(
					tempCSVfolder.getAbsolutePath() + File.separator + myProt.getFileName() + "_Chain_" + chain
							.getChainID() + "_zscoreCorrect.csv");

			File correctPositions = new File(
					tempCSVfolder.getAbsolutePath() + File.separator + myProt.getFileName() + "_Chain_" + chain
							.getChainID() +
							"_originalPositions.csv");

			writeMatrixToCSV(resultCSV, chain.intensityValueMatrix);
			writeMatrixToCSV(zscoreCSV, chain.zValueMatrix);
			writeTrueValueCSVs(zscoreCorrect, chain);
			writeMatrixToCSV(correctPositions, chain.originalPositions);

		}


	}

	private void writeTrueValueCSVs(File outputCSV, SimpleProtein.ProtChain chain) throws
			IOException {
		FileWriter FW = new FileWriter(outputCSV);

		String row = "";
		for (int j = 0; j < chain.originalAcidsScore.length; j++) {
			row += chain.originalPositions[j] + ", " + chain.originalAcidsScore[j] + "\n";
		}
		FW.write(row);

		FW.close();
	}


	private void writeMatrixToCSV(File outputCSV, double[][] matrix) throws IOException {
		FileWriter FW = new FileWriter(outputCSV);
		for (int i = 0; i < matrix[0].length; i++) {
			String row = "";
			for (int j = 0; j < matrix.length; j++) {
				row += matrix[j][i] + ", ";
			}
			row += "\n";
			FW.write(row);
		}
		FW.close();
	}

	private void writeMatrixToCSV(File outputCSV, double[] matrix) throws IOException {
		FileWriter FW = new FileWriter(outputCSV);

		String row = "";
		for (int j = 0; j < matrix.length; j++) {
			row += matrix[j] + "\n";
		}
		FW.write(row);

		FW.close();
	}

	private void writeMatrixToCSV(File outputCSV, Integer[] matrix) throws IOException {
		FileWriter FW = new FileWriter(outputCSV);

		String row = "";
		for (int j = 0; j < matrix.length; j++) {
			row += matrix[j] + "\n";
		}
		FW.write(row);

		FW.close();
	}

	private Integer[] getOriginalPositions(SimpleProtein myProt) throws InvalidPropertiesFormatException {
		List<Integer> positionArray = new ArrayList<>();
		for (SimpleProtein.ProtChain chain : myProt) {
			for (AminoAcid acid : chain) {
				positionArray.add(acidToIndex(acid.getName()));
			}
		}
		return positionArray.toArray(new Integer[positionArray.size()]);
	}


	public SimpleProtein getMyProt() {
		return myProt;
	}
}
