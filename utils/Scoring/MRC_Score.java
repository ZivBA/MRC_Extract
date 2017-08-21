package utils.Scoring;

import utils.ExtractMaxValue;
import utils.ScoreUtilities.MRC_Map_New;
import utils.ScoreUtilities.ScoringGeneralHelpers;
import utils.UtilExceptions.MissingChainID;
import utils.molecularElements.AminoAcid;
import utils.molecularElements.ProteinActions;
import utils.molecularElements.SimpleAtom;
import utils.molecularElements.SimpleProtein;
import utils.scwrlIntegration.SCWRLactions;

import static utils.ScoreUtilities.ScoringGeneralHelpers.csvToMatrix;
import static utils.ScoreUtilities.ScoringGeneralHelpers.debug;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import static org.apache.commons.lang3.StringUtils.substringBetween;
import static utils.molecularElements.ProteinActions.acidToIndex;

/**
 * Created by zivben on 09/08/15.
 */
public class MRC_Score {
	private List<Integer[]> originalPos;
	private MRC_Map_New myMap;
	private SimpleProtein myProt;
	public char requestedChain;
	private double[][] proteinIntensityValueMatrix;
	public ArrayList<String> logFile = new ArrayList<>();
	public ArrayList<File> toDelete = new ArrayList<>();
	private DateFormat df = new SimpleDateFormat("dd/MM/yy HH:mm:ss");


	public MRC_Score(MRC_Map_New myMap, SimpleProtein myProt) {
		this.myMap = myMap;
		this.myProt = myProt;
		proteinIntensityValueMatrix = new double[20][myProt.getNumChains()];
		
//		exMax();

	}

	public MRC_Score(String mapPath, String protPath) throws IOException {
//		System.out.println("Read PDB");
		this.myProt = new SimpleProtein(new File(protPath));
		if (!checkExistingCSVs()){
			this.myMap = new MRC_Map_New(mapPath);
			File pdbFile = new File(protPath);
			SCWRLactions.scwrlRunOnce(pdbFile,pdbFile);
			this.myProt = new SimpleProtein(new File(protPath));
		}
//		exMax();

		proteinIntensityValueMatrix = new double[20][myProt.getNumChains()];


	}

	public MRC_Score(String mapPath, String protPath, String requestedChain) throws IOException {
//		System.out.println("Read PDB");
		this.myProt = new SimpleProtein(new File(protPath));
		if (!checkExistingCSVs()){
			this.myMap = new MRC_Map_New(mapPath);
			File pdbFile = new File(protPath);
			SCWRLactions.scwrlRunOnce(pdbFile,pdbFile);
			this.myProt = new SimpleProtein(new File(protPath));
		}
		this.requestedChain = requestedChain.charAt(0);
		proteinIntensityValueMatrix = new double[20][myProt.getNumChains()];

	}

	public MRC_Score(SimpleProtein myProt, String requestedChain, MRC_Map_New myMap) throws IOException {
		this.myProt = myProt;
		this.myMap = myMap;
		
		this.requestedChain = requestedChain.trim().charAt(0);
		proteinIntensityValueMatrix = new double[20][myProt.getNumChains()];
		logFile.add("Started at: "+ df.format(new Date()));
	}

	private void exMax() {
		float[] maxValResult = ExtractMaxValue.getMaxValue(myMap);
		System.out.println(Arrays.toString(maxValResult));
		try {
			ExtractMaxValue.writeMarkerFile(myProt.getSource().getParent(), maxValResult);
		} catch (IOException e) {
			e.printStackTrace();
		}
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

		// if cannot read existing results, generate new files.
		if (!checkExistingCSVs()) {

			logFile.add("Starting protein scoring, saving original positions.");
			
			originalPos = myProt.getOriginalPositions();
			SimpleProtein backup = myProt;

			logFile.add("Stripping all amino acid resiudes and setting to ALA.");
			ProteinActions.stripAndAllALAToObject(myProt);
			logFile.add("Iterating all acid permutations and creating SCWRL input files");
			File scwrlOutput = ProteinActions.iterateAcids(myProt,requestedChain);

			File chainSubFolders[] = scwrlOutput.listFiles();
			for (File chain : chainSubFolders) {
				if (chain.isDirectory()) {
					if (chain.getAbsolutePath().endsWith(File.separator + requestedChain) || requestedChain == '\0') {
						SCWRLactions.genSCWRLforFolder(chain);
						if (!debug) {
							toDelete.add(chain);
						}
					}
				}
			}
			myProt = backup;
			processSCWRLfolder(scwrlOutput, requestedChain);


		}
		return proteinIntensityValueMatrix;
	}

	private boolean checkExistingCSVs() {

		//try reading CSVs
		try {
			File tempCSVfolder = ScoringGeneralHelpers.makeFolder(new File(myProt.getSource().getParent() + File
					.separator + "tempCSVs"));


			for (SimpleProtein.ProtChain chain : myProt) {
				String upperCaseFileName = myProt.getFileName().toUpperCase();
				if (chain.getChainID() == requestedChain || requestedChain == '\0') {
					File resultCSV = new File(tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
							+ "_" + chain.getChainID() + "_resultMatrix.csv");

					File backBoneIntMatrix = new File(
							tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
									+ "_" + chain.getChainID() + "_BBresultMatrix.csv");

					File correctPositions = new File(
							tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
									+ "_" + chain.getChainID() + "_originalPositions.csv");

					double[][] temp;
					temp = ScoringGeneralHelpers.csvToMatrix(resultCSV);
					if (temp != null) {
						chain.resIntensityValueMatrix = temp;
					} else {
						return false;
					}
					temp = csvToMatrix(backBoneIntMatrix);
					if (temp != null) {
						chain.backBoneIntensityValueMatrix = temp;
					} else {
						return false;
					}
					double[][] tempPositions = csvToMatrix(correctPositions);
					if (tempPositions == null)
						return false;
					chain.originalPositions = new Integer[tempPositions[0].length];
					for (int i = 0; i < tempPositions.length; i++) {
						for (int j = 0; j < tempPositions[i].length; j++) {
							chain.originalPositions[j] = (int) tempPositions[i][j];
						}
					}
				}
//				System.out.println("all good with existing CSVs");
			
			}
			return true;
		} catch (Exception e) {
			System.out.println("Existing CSVs not found or malformed, processing from scratch.");
			return false;
		}
		
	}

	private void processSCWRLfolder(File processingFolder, char requestedChain) throws IOException, MissingChainID {


		SimpleProtein tempProt;

		/**
		 * get input File object, read and convert to List of Files (PDB Iterations) for processing.
		 */
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
			logFile.add(ex.getLocalizedMessage());
			logFile.add(ex.getMessage());
			throw ex;
		}

		/**
		 * read each chain File object
		 */
		for (File chain : chainFolders) {
			try (DirectoryStream<Path> directoryStream = Files.newDirectoryStream(chain.toPath())) {
				for (Path path : directoryStream) {
					tempProt = new SimpleProtein(path.toFile());
					scoreSingleScwrl(tempProt);
				}

			} catch (IOException | MissingChainID ex) {
				logFile.add("array index OOB exception for file: " + chain.getAbsolutePath());
				logFile.add(ex.getMessage());
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
			logFile.add(ex.getLocalizedMessage());
			logFile.add(ex.getMessage());
			throw ex;
		}

		for (File chain : chainFolders) {
			try (DirectoryStream<Path> directoryStream = Files.newDirectoryStream(chain.toPath())) {
				for (Path path : directoryStream) {
					tempProt = new SimpleProtein(path.toFile());
					scoreSingleScwrl(tempProt);
				}
			} catch (ArrayIndexOutOfBoundsException ex) {
				System.out.println("array index OOB exception for file: " + chain.getAbsolutePath());
			} catch (IOException | MissingChainID ex) {
				logFile.add("array index OOB exception for file: " + chain.getAbsolutePath());
				logFile.add(ex.getMessage());
				throw ex;
			}


		}
	}

	private void scoreSingleScwrl(SimpleProtein tempProt) throws InvalidPropertiesFormatException, MissingChainID {
		int testResiduePosition = 0;
		try {
			 testResiduePosition = Integer.valueOf(substringBetween(tempProt.getFileName(), "_res_", "_to"));
		} catch (NumberFormatException e){
			System.err.println("Problem with some SCWRL file?");
			System.out.println("Problematic file is:");
			System.out.println(tempProt.getSource().getAbsolutePath());
			System.out.println("also this is the filename analized: "+tempProt.getFileName());
			throw e;
		}
		int testResidueIndex = acidToIndex(
				substringBetween(tempProt.getFileName(), "_to_", "_SCWRLed"));

		for (SimpleProtein.ProtChain chain : tempProt) {
			if (chain.getChainID() == requestedChain || requestedChain == '\0') {
				for (AminoAcid residue : chain) {

					if (residue.getPosition() == testResiduePosition && acidToIndex(residue.getName()) == testResidueIndex) {

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
						try {
							originalChain.resIntensityValueMatrix[acidToIndex(residue.getName())][residue.getPosition()] = resSum;
							originalChain.backBoneIntensityValueMatrix[acidToIndex(residue.getName())][residue.getPosition()] = backBoneSum;
						}
						catch (ArrayIndexOutOfBoundsException e){
							System.out.println("OOB exception at residue: "+ residue.getName() + " in position: "+ residue.getPosition());
						}
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
			logFile.add("Warning, PDB contains coordinate value outside MRC map scope.");
			logFile.add("Atom " + atom.getName() + " from residue " + atom.getaAcidName() + " at coordinates " +
					Arrays.toString(atom.getAtomCoords()) + " is the problem.");
			logFile.add("problematic line from PDB is: \n" + atom.getOriginalString() + "\n");
			atom.setAtomScore(0.0);
		}
		return atom.getAtomScore();
	}

	public void calcZvalue() throws InvalidPropertiesFormatException {

		for (SimpleProtein.ProtChain chain : myProt) {

			if (chain.getChainID() == requestedChain) {
				logFile.add("Calculating Z-Values for chain: " + chain.getChainID());
				calcMedian(chain);
				zValueHelper(chain);
				//			zValueHelperWithBackBone(chain);
			}
		}


	}

	private void calcMedian(SimpleProtein.ProtChain chain) {
//		List<List<Double>> falseValuesList = new LinkedList<>();
//		List<List<Double>> trueValuesList = new LinkedList<>();
		List<List<Double>> resValueList = new LinkedList<>();
		List<List<Double>> BBValueList = new LinkedList<>();

		for (int i = 0; i < chain.resIntensityValueMatrix.length; i++) {
//			falseValuesList.add(new LinkedList<Double>());
//			trueValuesList.add(new LinkedList<Double>());
			BBValueList.add(new LinkedList<Double>());
			resValueList.add(new LinkedList<Double>());
			for (int j = 0; j < chain.resIntensityValueMatrix[i].length; j++) {
//				double tmpScore = ;
//				if (chain.originalPositions[j] != i) {
//					falseValuesList.get(i).add(tmpScore);
//				} else {
//					trueValuesList.get(i).add(tmpScore);
//				}
				BBValueList.get(i).add(chain.backBoneIntensityValueMatrix[i][j]);
				resValueList.get(i).add(chain.resIntensityValueMatrix[i][j]);
				
			}
//			Collections.sort(falseValuesList.get(i));
//			Collections.sort(trueValuesList.get(i));
			Collections.sort(BBValueList.get(i));
			Collections.sort(resValueList.get(i));
		}
//		double[] falseValuesMedian = new double[falseValuesList.size()];
//		double[] trueValuesMedian = new double[trueValuesList.size()];
		double[] BBValuesMedian = new double[BBValueList.size()];
		double[] resValuesMedian = new double[resValueList.size()];

//
//		for (int i = 0; i < falseValuesList.size(); i++) {
//			if (falseValuesList.get(i).size() == 0) {
//				falseValuesMedian[i] = 0;
//			} else {
//				int middle = falseValuesList.get(i).size() / 2;
//				if (falseValuesList.get(i).size() % 2 == 1) {
//					falseValuesMedian[i] = falseValuesList.get(i).get(middle);
//				} else {
//					falseValuesMedian[i] = (falseValuesList.get(i).get(middle - 1) + falseValuesList.get(i).get(middle)) / 2.0;
//				}
//			}
//		}

	
		for (int i = 0; i < BBValueList.size(); i++) {
			if (BBValueList.get(i).size() == 0) {
				BBValuesMedian[i] = 0;
			} else {
				int middle = BBValueList.get(i).size() / 2;
				if (BBValueList.get(i).size() % 2 == 1) {
					BBValuesMedian[i] = BBValueList.get(i).get(middle);
				} else {
					BBValuesMedian[i] = (BBValueList.get(i).get(middle - 1) + BBValueList.get(i).get(middle)) / 2.0;
				}
			}
		}

		for (int i = 0; i < resValueList.size(); i++) {
			if (resValueList.get(i).size() == 0) {
				resValuesMedian[i] = 0;
			} else {
				int middle = resValueList.get(i).size() / 2;
				if (resValueList.get(i).size() % 2 == 1) {
					resValuesMedian[i] = resValueList.get(i).get(middle);
				} else {
					resValuesMedian[i] = (resValueList.get(i).get(middle - 1) + resValueList.get(i).get(middle)) / 2.0;
				}
			}
		}

//		chain.medianTrue = trueValuesMedian;
//		chain.medianFalse = falseValuesMedian;
		chain.allMedian = resValuesMedian;
		chain.backBoneMedian = BBValuesMedian;

//		for (int i = 0; i < chain.signalMaybe.length; i++) {
//			chain.signalMaybe[i] = chain.medianTrue[i] - chain.medianFalse[i];
//		}

	}
	/**
	 * calc zValue per position/AA combination. only for residue atoms, no BB
	 *
	 * @param chain
	 */
	private void zValueHelper(SimpleProtein.ProtChain chain) {
		double tempAvg[] = new double[20];
		double tempStD[] = new double[20];
		chain.allZvalueMatrix = new double[chain.resIntensityValueMatrix.length][chain.resIntensityValueMatrix[0]
				.length];
		chain.backBoneZvalueMatrix = new double[chain.resIntensityValueMatrix.length][chain.resIntensityValueMatrix[0]
				.length];
//
		//*************************************
		// redo all for backbone, just in case.
		//*************************************


		for (int i = 0; i < chain.backBoneIntensityValueMatrix.length; i++) {
			for (int j = 0; j < chain.backBoneIntensityValueMatrix[i].length; j++) {
				tempAvg[i] += chain.backBoneIntensityValueMatrix[i][j];
			}
			tempAvg[i] = tempAvg[i] / chain.backBoneIntensityValueMatrix[i].length;

			//calc standard deviation for each column
			for (int j = 0; j < chain.backBoneIntensityValueMatrix[i].length; j++) {
				tempStD[i] += ((chain.backBoneIntensityValueMatrix[i][j] - tempAvg[i]) * (chain.backBoneIntensityValueMatrix[i][j] - tempAvg[i]));
			}

			tempStD[i] = (double) Math.round((Math.sqrt(tempStD[i] / chain.backBoneIntensityValueMatrix[i].length)) * 10000000d) / 10000000d;
		}

		// calc Z-Value for each discrete acid in every position
		for (int i = 0; i < chain.backBoneIntensityValueMatrix.length; i++) {
			for (int j = 0; j < chain.backBoneIntensityValueMatrix[i].length; j++) {

				Double tmpScore = (double) Math.round(
						(((chain.backBoneIntensityValueMatrix[i][j]) - chain.backBoneMedian[i]) / tempStD[i]) * 10000000d) / 10000000d;
				if (tempStD[i] < (0.000009d)) {
					chain.backBoneZvalueMatrix[i][j] = 0.0d;
				} else {

					if (tmpScore.isNaN()) {
						tmpScore = 0.0;
					}
					// add score to original acid score only if right position
					if (chain.originalPositions[j] == i) {
						chain.BBTrueZvalue[j] = tmpScore > 4 ? 4 : (tmpScore < -4 ? -4 : tmpScore);;
					}
					// add all scores to the allZvalueMatrix.
					chain.backBoneZvalueMatrix[i][j]  = tmpScore > 4 ? 4 : (tmpScore < -4 ? -4 : tmpScore);

				}
			}
		}

		// calc tempAvg per column in the intensity value matrix ( avarage of scores per amino acid in every pos)

		tempAvg = new double[20];
		tempStD = new double[20];

		// calc tempAvg per column in the intensity value matrix ( avarage of scores per amino acid in every pos)
		for (int i = 0; i < chain.resIntensityValueMatrix.length; i++) {
			for (int j = 0; j < chain.resIntensityValueMatrix[i].length; j++) {
				tempAvg[i] += chain.resIntensityValueMatrix[i][j];
			}
			tempAvg[i] = tempAvg[i] / chain.resIntensityValueMatrix[i].length;

			//calc standard deviation for each column
			for (int j = 0; j < chain.resIntensityValueMatrix[i].length; j++) {
				double tmpResScore = chain.resIntensityValueMatrix[i][j];
				tempStD[i] += ((tmpResScore - tempAvg[i]) * (tmpResScore - tempAvg[i]));
			}
			tempStD[i] = (double) Math.round((Math.sqrt(tempStD[i] / chain.resIntensityValueMatrix[i].length)) * 10000000d) / 10000000d;
		}

		// calc Z-Value for each discrete acid in every position
		for (int i = 0; i < chain.resIntensityValueMatrix.length; i++) {
			for (int j = 0; j < chain.resIntensityValueMatrix[i].length; j++) {
				double tmpResScore = chain.resIntensityValueMatrix[i][j];
				if (tempStD[i] < (0.000009d)) {
					chain.allZvalueMatrix[i][j] = 0.0d;
				} else {
					Double tmpScore = (double) Math.round(
							((tmpResScore - chain.allMedian[i]) / tempStD[i]) * 10000000d) / 10000000d;
					if (tmpScore.isNaN()) {
						tmpScore = 0.0;
					}
					// add score to original acid score only if right position
					if (chain.originalPositions[j] == i) {
						chain.trueZvalues[j] = tmpScore > 4 ? 4 : (tmpScore < -4 ? -4 : tmpScore);;
					}
					chain.allZvalueMatrix[i][j] = tmpScore > 4 ? 4 : (tmpScore < -4 ? -4 : tmpScore);;

				}
			}
		}



	}

//	/**
//	 * calc zValue with backbone weighting or addition (temp method).
//	 *
//	 * @param chain
//	 */

//	private void zValueHelperWithBackBone(SimpleProtein.ProtChain chain) {
//		double tempAvg[] = new double[20];
//		double tempStD[] = new double[20];
//		chain.allZvalueMatrix = new double[chain.resIntensityValueMatrix.length][chain.resIntensityValueMatrix[0]
//				.length];
//		chain.backBoneZvalueMatrix = new double[chain.resIntensityValueMatrix.length][chain.resIntensityValueMatrix[0]
//				.length];
//
//		//*************************************
//		//       for backbone atoms.
//		//*************************************
//
//		/*
//		calc avarage and std deviation per column (AA type).
//		 */
//		for (int i = 0; i < chain.backBoneIntensityValueMatrix.length; i++) {
//			for (int j = 0; j < chain.backBoneIntensityValueMatrix[i].length; j++) {
//				tempAvg[i] += chain.backBoneIntensityValueMatrix[i][j];
//			}
//			tempAvg[i] = tempAvg[i] / chain.backBoneIntensityValueMatrix[i].length;
//
//			//calc standard deviation for each column
//			for (int j = 0; j < chain.backBoneIntensityValueMatrix[i].length; j++) {
//				tempStD[i] += Math.pow(chain.backBoneIntensityValueMatrix[i][j] - tempAvg[i], 2);
//			}
//			tempStD[i] = Math.sqrt(tempStD[i] / chain.backBoneIntensityValueMatrix[i].length);
//		}
//
//		// calc Z-Value for each discrete acid in every position
//		for (int i = 0; i < chain.backBoneIntensityValueMatrix.length; i++) {
//			for (int j = 0; j < chain.backBoneIntensityValueMatrix[i].length; j++) {
//
//				Double tmpScore = (chain.backBoneIntensityValueMatrix[i][j] - tempAvg[i]) / tempStD[i];
//				if (tmpScore.isNaN()) {
//					tmpScore = 0.0;
//				}
//				// add score to original acid score only if right position
//				if (chain.originalPositions[j] == i) {
//					chain.BBTrueZvalue[j] = tmpScore;
//				}
//				// add all scores to the allZvalueMatrix.
//				chain.backBoneZvalueMatrix[i][j] = tmpScore;
//
//			}
//		}
//
//
//		//*************************************
//		//       for residue atoms.
//		//*************************************
//
//		tempAvg = new double[20];
//		tempStD = new double[20];
//		for (int i = 0; i < chain.resIntensityValueMatrix.length; i++) {
//			for (int j = 0; j < chain.resIntensityValueMatrix[i].length; j++) {
//
//
//				// weight the intensity value by multiplying with BB intensity value.
//				//				chain.resIntensityValueMatrix[i][j] *= chain.backBoneIntensityValueMatrix[i][j];
//
//				tempAvg[i] += (chain.resIntensityValueMatrix[i][j]);
//			}
//			tempAvg[i] = tempAvg[i] / chain.resIntensityValueMatrix[i].length;
//
//			//calc standard deviation for each column
//			for (int j = 0; j < chain.resIntensityValueMatrix[i].length; j++) {
//				tempStD[i] += Math.pow((chain.resIntensityValueMatrix[i][j]) - tempAvg[i], 2);
//			}
//			tempStD[i] = Math.sqrt(tempStD[i] / chain.resIntensityValueMatrix[i].length);
//		}
//
//		// calc Z-Value for each discrete acid in every position
//		for (int i = 0; i < chain.resIntensityValueMatrix.length; i++) {
//			for (int j = 0; j < chain.resIntensityValueMatrix[i].length; j++) {
//
//				Double tmpScore = ((chain.resIntensityValueMatrix[i][j]) -
//						tempAvg[i]) / tempStD[i];
//				if (tmpScore.isNaN()) {
//					tmpScore = 0.0;
//				}
//				// add score to original acid score only if right position
//				if (chain.originalPositions[j] == i) {
//					chain.trueZvalues[j] = tmpScore;
//				}
//				chain.allZvalueMatrix[i][j] = tmpScore;
//
//			}
//		}
//
//
//	}

	public void createCSVs() throws IOException {

		logFile.add("Creating CSV Files");
		File tempCSVfolder = ScoringGeneralHelpers.makeFolder(new File(myProt.getSource().getParent() + File.separator + "tempCSVs"));
		String upperCaseFileName = myProt.getFileName().toUpperCase();
		SimpleProtein.ProtChain chain = myProt.getChain(requestedChain);


		// run zvaluematrix through de-negativeation vector (try to make all values positive)


		File resultCSV = new File(tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
				+ "_" + chain.getChainID() + "_resultMatrix.csv");

		File backBoneIntMatrix = new File(tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
				+ "_" + chain.getChainID() + "_BBresultMatrix.csv");

		File zscoreCSV = new File(tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
				+ "_" + chain.getChainID() + "_allZscore.csv");

		File zscoreBBCSV = new File(tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
				+ "_" + chain.getChainID() + "_allZscoreWtBB.csv");
		File zscoreOnlyFalseCSV = new File(tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
				+ "_" + chain.getChainID() + "_onlyFalseZscore.csv");

		File trueMedian = new File(tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
				+ "_" + chain.getChainID() + "_trueMedian.csv");

		File falseMedian = new File(tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
				+ "_" + chain.getChainID() + "_falseMedian.csv");

		File signalMaybe = new File(tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
				+ "_" + chain.getChainID() + "_signalMaybe.csv");

		File combinedMatrixNoVec = new File(tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
				+ "_" + chain.getChainID() + "_profileNoVec.txt");
		
		
		File combinedMatrixLatestVec = new File(tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
				+ "_" + chain.getChainID() + "_profileLatestVec.txt");

//		File combinedMatrixLatestVecWeightBB2 = new File(tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
//				+ "_" + chain.getChainID() + "_profileLatestVec_weightedBB_2.txt");
//
//		File combinedMatrixLatestVecWeightBB5 = new File(tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
//				+ "_" + chain.getChainID() + "_profileLatestVec_weightedBB_5.txt");
//		File combinedMatrixLatestVecWeightBB10 = new File(tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName
//				+ "_" + chain.getChainID() + "_profileLatestVec_weightedBB_10.txt");

		File zscoreCorrect = new File(
				tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName + "_" + chain
						.getChainID() + "_zscoreCorrect.csv");

		File correctPositions = new File(
				tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName + "_" + chain
						.getChainID() +
						"_originalPositions.csv");

		File backBoneZscore = new File(
				tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName + "_" + chain
						.getChainID() + "_backboneZscore.csv");


		File logFileTarget = new File(
				tempCSVfolder.getAbsolutePath() + File.separator + upperCaseFileName + "_" +
						chain.getChainID() + "_log.txt");

		writeMatrixToCSV(trueMedian, chain.medianTrue);
		writeMatrixToCSV(falseMedian, chain.medianFalse);
		writeMatrixToCSV(signalMaybe, chain.signalMaybe);

		// generate profile files per vector
		//			writeNewMatrixFormat(combinedMatrixVec4, chain, ScoringGeneralHelpers.vector4);
		//			writeNewMatrixFormat(combinedMatrixVec3, chain, ScoringGeneralHelpers.vector3);
		writeNewMatrixFormat(combinedMatrixLatestVec, chain, ScoringGeneralHelpers.latestVector, 0);
		writeNewMatrixFormat(combinedMatrixLatestVec, chain, ScoringGeneralHelpers.latestVector, 12);
		writeNewMatrixFormat(combinedMatrixLatestVec, chain, ScoringGeneralHelpers.latestVector, 8);
		writeNewMatrixFormat(combinedMatrixLatestVec, chain, ScoringGeneralHelpers.latestVector, 4);
		writeNewMatrixFormat(combinedMatrixLatestVec, chain, ScoringGeneralHelpers.latestVector, 2);
		writeNewMatrixFormat(combinedMatrixLatestVec, chain, ScoringGeneralHelpers.latestVector, 1);
		writeNewMatrixFormat(combinedMatrixNoVec, chain, ScoringGeneralHelpers.normalVector, 0);
		writeNewMatrixFormat(combinedMatrixNoVec, chain, ScoringGeneralHelpers.normalVector, 12);
		writeNewMatrixFormat(combinedMatrixNoVec, chain, ScoringGeneralHelpers.normalVector, 8);
		writeNewMatrixFormat(combinedMatrixNoVec, chain, ScoringGeneralHelpers.normalVector, 4);
		writeNewMatrixFormat(combinedMatrixNoVec, chain, ScoringGeneralHelpers.normalVector, 2);
		writeNewMatrixFormat(combinedMatrixNoVec, chain, ScoringGeneralHelpers.normalVector, 1);

		writeMatrixToCSV(resultCSV, chain.resIntensityValueMatrix);
		writeMatrixToCSV(backBoneIntMatrix, chain.backBoneIntensityValueMatrix);
		writeFalseZvalueResults(zscoreOnlyFalseCSV, chain);
		writeMatrixToCSV(zscoreCSV, chain.allZvalueMatrix);
		writeMatrixToCSV(backBoneZscore, chain.BBTrueZvalue);
		writeMatrixToCSV(correctPositions, chain.originalPositions);
		writeTrueValueCSVs(zscoreCorrect, chain);
		writeBBzScoreToCSV(zscoreBBCSV, chain, 5);

		logFile.add("Finished at: " + df.format(new Date()));
		wirteLogToTxtFile(logFileTarget, logFile);




	}

	private void wirteLogToTxtFile(File logFileTarget, ArrayList<String> logFile) throws IOException {
		FileWriter FW = new FileWriter(logFileTarget);
		for (String line : logFile) {
			FW.write(line + "\n");
		}
		FW.close();
	}

	private void writeNewMatrixFormat(File outputCSV, SimpleProtein.ProtChain chain, double[] vector, double weight)
			throws
			IOException {
		
//		calcZvalue();

		double[][] zvaltemp = new double[chain.allZvalueMatrix.length][chain.allZvalueMatrix[0].length];
		for (int i = 0; i < zvaltemp.length; i++) {
			for (int j = 0; j < zvaltemp[0].length; j++) {
				zvaltemp[i][j] = chain.allZvalueMatrix[i][j];
			}
		}
		


		ScoringGeneralHelpers.multiplyMatrixByVector(zvaltemp, vector);
		
		FileWriter FW = new FileWriter(outputCSV);
		for (int i = 0; i < zvaltemp[0].length; i++) {
			String row = "";
			row += chain.getAcidSequenceID(i) + "\t";
			row += chain.originalPositions[i] + "\t";
			for (int j = 0; j < chain.allZvalueMatrix.length; j++) {
				row += (zvaltemp[j][i]) + "\t";
			}
			
			row += "\n";
			FW.write(row);
		}
		FW.close();
		
		zvaltemp = new double[chain.allZvalueMatrix.length][chain.allZvalueMatrix[0].length];
		for (int i = 0; i < zvaltemp.length; i++) {
			for (int j = 0; j < zvaltemp[0].length; j++) {
				zvaltemp[i][j] = chain.allZvalueMatrix[i][j]*(1+(weight*(chain.backBoneZvalueMatrix[i][j] / 4)));
			}
		}
		
		ScoringGeneralHelpers.multiplyMatrixByVector(zvaltemp, vector);
		
		FileWriter FWWBB = new FileWriter(new File(outputCSV.getAbsolutePath().replace(".txt", "_weightedBB_" + weight + ".txt")));

		for (int i = 0; i < zvaltemp[0].length; i++) {
			String row = "";
			row += chain.getAcidSequenceID(i) + "\t";
			row += chain.originalPositions[i] + "\t";
			for (int j = 0; j < zvaltemp.length; j++) {
				row += (zvaltemp[j][i]) + "\t";
			}

			row += "\n";
			FWWBB.write(row);
		}
		FWWBB.close();


	}

	private void writeFalseZvalueResults(File outputCSV, SimpleProtein.ProtChain chain) throws IOException {
		FileWriter FW = new FileWriter(outputCSV);
		for (int i = 0; i < chain.allZvalueMatrix[0].length; i++) {
			String row = "";
			for (int j = 0; j < chain.allZvalueMatrix.length; j++) {
				if (chain.originalPositions[i] == j) {
					row += 0.0;
				} else {
					row += chain.allZvalueMatrix[j][i];
				}
				if (j != chain.allZvalueMatrix.length - 1) {
					row += ", ";
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

			//			row += chain.BBTrueZvalue[j] + "\n";
		}
		FW.write(row);

		FW.close();
	}


	private void writeMatrixToCSV(File outputCSV, double[][] matrix) throws IOException {
		FileWriter FW = new FileWriter(outputCSV);
		for (int i = 0; i < matrix[0].length; i++) {
			String row = "";
			for (int j = 0; j < matrix.length; j++) {
				row += matrix[j][i];
				if (j != matrix.length - 1) {
					row += ", ";
				}
			}
			row += "\n";
			FW.write(row);
		}
		FW.close();
	}

	private void writeBBzScoreToCSV(File outputCSV, SimpleProtein.ProtChain chain, int weight) throws IOException {


		FileWriter FW = new FileWriter(outputCSV);
		for (int i = 0; i < chain.allZvalueMatrix[0].length; i++) {
			String row = "";
			for (int j = 0; j < chain.allZvalueMatrix.length; j++) {
				row += (chain.allZvalueMatrix[j][i] * (1 + (chain.backBoneZvalueMatrix[j][i] / weight)));
				if (j != chain.allZvalueMatrix.length - 1) {
					row += ", ";
				}
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


	public SimpleProtein getMyProt() {
		return myProt;
	}

	public int[] getAcidDist() throws InvalidPropertiesFormatException {
		return (myProt.acidDist != null ? myProt.acidDist : myProt.calcAcidDist());
	}
}
