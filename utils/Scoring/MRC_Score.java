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
import java.util.*;

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
						originalChain.resIntensityValueMatrix[acidToIndex(residue.getName())][residue.getSeqNum()
								- originalChain.getSequenceBias()] = resSum;
						originalChain.backBoneIntensityValueMatrix[acidToIndex(residue.getName())][residue.getSeqNum
								() - originalChain.getSequenceBias()] = backBoneSum;

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
			calcMedian(chain);
		}


	}

	private void calcMedian(SimpleProtein.ProtChain chain) {
		List<List<Double>> falseValuesList = new LinkedList<>();
		List<List<Double>> trueValuesList = new LinkedList<>();

		for (int i = 0; i < chain.allZvalueMatrix.length; i++) {
			falseValuesList.add(new LinkedList<Double>());
			trueValuesList.add(new LinkedList<Double>());
			for (int j = 0; j < chain.allZvalueMatrix[i].length; j++) {
				if (chain.originalPositions[j] != i) {
					falseValuesList.get(i).add(chain.allZvalueMatrix[i][j]);
				} else {
					trueValuesList.get(i).add(chain.allZvalueMatrix[i][j]);
				}
			}
			Collections.sort(falseValuesList.get(i));
			Collections.sort(trueValuesList.get(i));
		}
		double[] falseValuesMedian = new double[falseValuesList.size()];
		double[] trueValuesMedian = new double[trueValuesList.size()];


		for (int i = 0; i < falseValuesList.size(); i++) {
			if (falseValuesList.get(i).size() == 0) {
				falseValuesMedian[i] = 0;
			} else {
				int middle = falseValuesList.get(i).size() / 2;
				if (falseValuesList.get(i).size() % 2 == 1) {
					falseValuesMedian[i] = falseValuesList.get(i).get(middle);
				} else {
					falseValuesMedian[i] = (falseValuesList.get(i).get(middle - 1) + falseValuesList.get(i).get(
							middle)) / 2.0;
				}
			}
		}


		for (int i = 0; i < trueValuesList.size(); i++) {
			if (trueValuesList.get(i).size() == 0) {
				trueValuesMedian[i] = 0;
			} else {
				int middle = trueValuesList.get(i).size() / 2;
				if (trueValuesList.get(i).size() % 2 == 1) {
					trueValuesMedian[i] = trueValuesList.get(i).get(middle);
				} else {
					trueValuesMedian[i] = (trueValuesList.get(i).get(middle - 1) + trueValuesList.get(i).get(
							middle)) / 2.0;
				}
			}
		}

		chain.medianTrue = trueValuesMedian;
		chain.medianFalse = falseValuesMedian;

		for (int i = 0; i < chain.signalMaybe.length; i++) {
			chain.signalMaybe[i] = chain.medianTrue[i] - chain.medianFalse[i];
		}



		/*double[][] falseValuesCopy = Arrays.copyOf(chain.allZvalueMatrix, chain.allZvalueMatrix.length);
		for (int i = 0; i < falseValuesCopy.length; i++) {
			Arrays.sort(falseValuesCopy[i]);

			int middle = falseValuesCopy[i].length / 2;
			if (falseValuesCopy[i].length % 2 == 1) {
				chain.medianFalse[i] = falseValuesCopy[i][middle];
			} else {
				chain.medianFalse[i] = (falseValuesCopy[i][middle - 1] + falseValuesCopy[i][middle]) / 2.0;
			}

		}
		double[] tempArray = Arrays.copyOf(chain.trueZvalues, chain.trueZvalues.length);
		Arrays.sort(tempArray);

		int middle = tempArray.length / 2;
		if (tempArray.length % 2 == 1) {
			chain.medianTrue = tempArray[middle];
		} else {
			chain.medianTrue = (tempArray[middle - 1] + tempArray[middle]) / 2.0;
		}*/


	}

	private void zValueHelper(SimpleProtein.ProtChain chain) {
		double tempAvg[] = new double[20];
		double tempStD[] = new double[20];
		chain.allZvalueMatrix = new double[chain.resIntensityValueMatrix.length][chain.resIntensityValueMatrix[0]
				.length];
		chain.backBoneZvalueMatrix = new double[chain.resIntensityValueMatrix.length][chain.resIntensityValueMatrix[0]
				.length];


		// calc tempAvg per column in the intensity value matrix ( avarage of scores per amino acid in every pos)
		for (int i = 0; i < chain.resIntensityValueMatrix.length; i++) {
			for (int j = 0; j < chain.resIntensityValueMatrix[i].length; j++) {
				tempAvg[i] += chain.resIntensityValueMatrix[i][j];
			}
			tempAvg[i] = tempAvg[i] / chain.resIntensityValueMatrix[i].length;

			//calc standard deviation for each column
			for (int j = 0; j < chain.resIntensityValueMatrix[i].length; j++) {
				tempStD[i] += Math.pow(chain.resIntensityValueMatrix[i][j] - tempAvg[i], 2);
			}
			tempStD[i] = Math.sqrt(tempStD[i] / chain.resIntensityValueMatrix[i].length);
		}

		// calc Z-Value for each discrete acid in every position
		for (int i = 0; i < chain.resIntensityValueMatrix.length; i++) {
			for (int j = 0; j < chain.resIntensityValueMatrix[i].length; j++) {

				double tmpScore = (chain.resIntensityValueMatrix[i][j] - tempAvg[i]) / tempStD[i];

				// add score to original acid score only if right position
				if (chain.originalPositions[j] == i) {
					chain.trueZvalues[j] = tmpScore;
				}
				chain.allZvalueMatrix[i][j] = tmpScore;

			}
		}

		//*************************************
		// redo all for backbone, just in case.
		//*************************************
		tempAvg = new double[20];
		tempStD = new double[20];

		for (int i = 0; i < chain.backBoneIntensityValueMatrix.length; i++) {
			for (int j = 0; j < chain.backBoneIntensityValueMatrix[i].length; j++) {
				tempAvg[i] += chain.backBoneIntensityValueMatrix[i][j];
			}
			tempAvg[i] = tempAvg[i] / chain.backBoneIntensityValueMatrix[i].length;

			//calc standard deviation for each column
			for (int j = 0; j < chain.backBoneIntensityValueMatrix[i].length; j++) {
				tempStD[i] += Math.pow(chain.backBoneIntensityValueMatrix[i][j] - tempAvg[i], 2);
			}
			tempStD[i] = Math.sqrt(tempStD[i] / chain.backBoneIntensityValueMatrix[i].length);
		}

		// calc Z-Value for each discrete acid in every position
		for (int i = 0; i < chain.backBoneIntensityValueMatrix.length; i++) {
			for (int j = 0; j < chain.backBoneIntensityValueMatrix[i].length; j++) {

				double tmpScore = (chain.backBoneIntensityValueMatrix[i][j] - tempAvg[i]) / tempStD[i];

				// add score to original acid score only if right position
				if (chain.originalPositions[j] == i) {
					chain.backBoneZvalue[j] = tmpScore;
				}
				// add all scores to the allZvalueMatrix.
				chain.backBoneZvalueMatrix[i][j] = tmpScore;

			}
		}


		//TODO redundant code?
		//		Integer[] suspectedCorrectPositions = chain.getOriginalPositions();
		//
		//		for (int i = 0; i < suspectedCorrectPositions.length; i++) {
		//			chain.zValueCorrect[i] = chain.allZvalueMatrix[suspectedCorrectPositions[i]][i];
		//		}


	}

	public void createCSVs() throws IOException {

		System.out.println("Creating CSV Files");
		File tempCSVfolder = ScoringGeneralHelpers.makeFolder(new File(myProt.getSource().getParent() + File
				.separator + "tempCSVs"));


		for (SimpleProtein.ProtChain chain : myProt) {

			// run zvaluematrix through de-negativeation vector (try to make all values positive)
			ScoringGeneralHelpers.multiplyMatrixByVector(chain.allZvalueMatrix);
			ScoringGeneralHelpers.multiplyMatrixByVector(chain);


			File resultCSV = new File(tempCSVfolder.getAbsolutePath() + File.separator + myProt.getFileName()
					+ "_" + chain.getChainID() + "_resultMatrix.csv");

			File zscoreCSV = new File(tempCSVfolder.getAbsolutePath() + File.separator + myProt.getFileName()
					+ "_" + chain.getChainID() + "_allZscore.csv");
			File zscoreOnlyFalseCSV = new File(tempCSVfolder.getAbsolutePath() + File.separator + myProt.getFileName()
					+ "_" + chain.getChainID() + "_onlyFalseZscore.csv");

			File trueMedian = new File(tempCSVfolder.getAbsolutePath() + File.separator + myProt.getFileName()
					+ "_" + chain.getChainID() + "_trueMedian.csv");

			File falseMedian = new File(tempCSVfolder.getAbsolutePath() + File.separator + myProt.getFileName()
					+ "_" + chain.getChainID() + "_falseMedian.csv");

			File signalMaybe = new File(tempCSVfolder.getAbsolutePath() + File.separator + myProt.getFileName()
					+ "_" + chain.getChainID() + "_signalMaybe.csv");

			File combinedMatrix = new File(tempCSVfolder.getAbsolutePath() + File.separator + myProt.getFileName()
					+ "_" + chain.getChainID() + "_allWithSequence.csv");




			File zscoreCorrect = new File(
					tempCSVfolder.getAbsolutePath() + File.separator + myProt.getFileName() + "_" + chain
							.getChainID() + "_zscoreCorrect.csv");

			File correctPositions = new File(
					tempCSVfolder.getAbsolutePath() + File.separator + myProt.getFileName() + "_" + chain
							.getChainID() +
							"_originalPositions.csv");

			File backBoneMatrix = new File(
					tempCSVfolder.getAbsolutePath() + File.separator + myProt.getFileName() + "_" + chain
							.getChainID() + "_backboneZscore.csv");

			writeMatrixToCSV(trueMedian, chain.medianTrue);
			writeMatrixToCSV(falseMedian, chain.medianFalse);
			writeMatrixToCSV(signalMaybe, chain.signalMaybe);

			writeNewMatrixFormat(combinedMatrix, chain);

			writeMatrixToCSV(resultCSV, chain.resIntensityValueMatrix);
			writeFalseZvalueResults(zscoreOnlyFalseCSV, chain);
			writeMatrixToCSV(zscoreCSV, chain.allZvalueMatrix);
			writeMatrixToCSV(backBoneMatrix, chain.backBoneZvalueMatrix);
			writeMatrixToCSV(correctPositions, chain.originalPositions);
			writeTrueValueCSVs(zscoreCorrect, chain);

		}


	}

	private void writeNewMatrixFormat(File outputCSV, SimpleProtein.ProtChain chain) throws IOException {


		FileWriter FW = new FileWriter(outputCSV);
		for (int i = 0; i < chain.allZvalueMatrix[0].length; i++) {
			String row = "";
			row += chain.getAcidSequenceID(i) + ", ";
			row += chain.originalPositions[i] + ", ";
			for (int j = 0; j < chain.allZvalueMatrix.length; j++) {
				row += chain.allZvalueMatrix[j][i] + ", ";
			}

			row += "\n";
			FW.write(row);
		}
		FW.close();


	}

	private void writeFalseZvalueResults(File outputCSV, SimpleProtein.ProtChain chain) throws IOException {
		FileWriter FW = new FileWriter(outputCSV);
		for (int i = 0; i < chain.allZvalueMatrix[0].length; i++) {
			String row = "";
			for (int j = 0; j < chain.allZvalueMatrix.length; j++) {
				if (chain.originalPositions[i] == j) {
					row += 0 + ", ";
				} else {
					row += chain.allZvalueMatrix[j][i] + ", ";
				}
			}
			row += "\n";
			FW.write(row);
		}
		FW.close();
	}

	private void writeTrueValueCSVs(File outputCSV, SimpleProtein.ProtChain chain) throws
			IOException {
		FileWriter FW = new FileWriter(outputCSV);
		String row;

		// uncomment for table headers.
		//		row = "Acid Sequence Id, Single Letter Name, Res zScore, Backbone zScore \n";
		//		FW.write(row);
		row = "";
		for (int j = 0; j < chain.trueZvalues.length; j++) {
			row += chain.getAcidSequenceID(j) + ", ";
			row += chain.originalPositions[j] + ", ";
			row += chain.trueZvalues[j] + "\n";

			//			row += chain.backBoneZvalue[j] + "\n";
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
