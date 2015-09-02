package utils.Scoring;

import utils.ScoreUtilities.MRC_Map_New;
import utils.ScoreUtilities.ScoringGeneralHelpers;
import utils.molecularElements.AminoAcid;
import utils.molecularElements.ProteinActions;
import utils.molecularElements.SimpleAtom;
import utils.molecularElements.SimpleProtein;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.InvalidPropertiesFormatException;
import java.util.List;

import static org.apache.commons.lang3.StringUtils.substringAfterLast;
import static org.apache.commons.lang3.StringUtils.substringBetween;
import static utils.molecularElements.ProteinActions.acidToIndex;

/**
 * Created by zivben on 09/08/15.
 */
public class MRC_Score {
	Integer[] originalPos;
	private MRC_Map_New myMap;
	private SimpleProtein myProt;
	private char requestedChain;
	private double[][] intensityValueMatrix;
	private double[][] zValueMatrix;
	private double[] zValueCorrect;
	private double[] originalAcidsScore;

	public MRC_Score(MRC_Map_New myMap, SimpleProtein myProt) {
		this.myMap = myMap;
		this.myProt = myProt;
		intensityValueMatrix = new double[20][myProt.getLegnth()];
		originalAcidsScore = new double[myProt.getLegnth()];

	}

	public MRC_Score(String mapPath, String protPath) throws IOException {
		this(new MRC_Map_New(mapPath), new SimpleProtein(new File(protPath)));
	}

	public MRC_Score(String mapPath, String protPath, String requestedChain) throws IOException {
		this(new MRC_Map_New(mapPath), new SimpleProtein(new File(protPath), requestedChain));
		this.requestedChain = requestedChain.charAt(0);

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
	public double[][] scoreProtein() throws IOException {
		System.out.println("Starting protein scoring, saving original positions.");
		myProt.saveOriginalPositions();
		originalPos = myProt.getOriginalPositions();

		System.out.println("Stripping all amino acid resiudes and setting to ALA.");
		ProteinActions.stripAndAllALAToObject(myProt);
		System.out.println("Iterating all acid permutations and creating SCWRL input files");
		File scwrlOutput = ProteinActions.iterateAndScwrl(myProt);

		processSCWRLfolder(scwrlOutput);

		return intensityValueMatrix;
	}

	private void processSCWRLfolder(File processingFolder) throws IOException {

		SimpleProtein tempProt;

		List<File> fileNames = new ArrayList<>();
		try (DirectoryStream<Path> directoryStream = Files.newDirectoryStream(processingFolder.toPath())) {
			for (Path path : directoryStream) {
				fileNames.add(path.toFile());
			}
		} catch (IOException ex) {
			throw ex;
		}

		for (File fileName : fileNames) {
			tempProt = new SimpleProtein(fileName);
			scoreSingleScwrl(tempProt);

		}
	}

	private void scoreSingleScwrl(SimpleProtein tempProt) throws InvalidPropertiesFormatException {


		int testResiduePosition = Integer.valueOf(substringBetween(tempProt.getFileName(), "_res_", "_"));
		int testResidueIndex = acidToIndex(
				substringAfterLast(tempProt.getFileName(), "_"));

		for (SimpleProtein.ProtChain chain : tempProt) {
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

					intensityValueMatrix[acidToIndex(
							residue.getName())][residue.getSeqNum() - myProt.getSequenceBias()] =
							/*backBoneSum  +*/ resSum;


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
			atom.setAtomScore(0.0);
		}
		return atom.getAtomScore();
	}

	public void calcZvalue() throws InvalidPropertiesFormatException {
		System.out.println("Calculating Z-Values");
		double tempAvg[] = new double[20];
		double tempStD[] = new double[20];
		zValueMatrix = new double[intensityValueMatrix.length][intensityValueMatrix[0].length];
		zValueCorrect = new double[myProt.getLegnth()];


		// calc tempAvg per column in the intensity value matrix ( avarage of scores per amino acid in every pos)
		for (int i = 0; i < intensityValueMatrix.length; i++) {
			for (int j = 0; j < intensityValueMatrix[i].length; j++) {
				tempAvg[i] += intensityValueMatrix[i][j];
			}
			tempAvg[i] = tempAvg[i] / intensityValueMatrix[i].length;

			//calc standard deviation for each column
			for (int j = 0; j < intensityValueMatrix[i].length; j++) {
				tempStD[i] += Math.pow(intensityValueMatrix[i][j] - tempAvg[i], 2);
			}
			tempStD[i] = Math.sqrt(tempStD[i] / intensityValueMatrix[i].length);
		}

		// calc Z-Value for each discrete acid in every position
		for (int i = 0; i < intensityValueMatrix.length; i++) {
			for (int j = 0; j < intensityValueMatrix[i].length; j++) {

				if (originalPos[j] == i) {
					originalAcidsScore[j] = (intensityValueMatrix[i][j] - tempAvg[i]) / tempStD[i];
					zValueMatrix[i][j] = 0;
				} else {
					zValueMatrix[i][j] = (intensityValueMatrix[i][j] - tempAvg[i]) / tempStD[i];
				}

			}
		}

		Integer[] suspectedCorrectPositions = myProt.getOriginalPositions();

		for (int i = 0; i < suspectedCorrectPositions.length; i++) {
			zValueCorrect[i] = zValueMatrix[suspectedCorrectPositions[i]][i];
		}

	}

	public void createCSVs() throws IOException {

		System.out.println("Creating CSV Files");
		File resultCSV = new File(myProt.getSource().getAbsolutePath().substring(0, myProt.getSource()
				.getAbsolutePath().indexOf(".pdb")) + "_resultMatrix.csv");
		File zscoreCSV = new File(myProt.getSource().getAbsolutePath().substring(0, myProt.getSource()
				.getAbsolutePath().indexOf(".pdb")) + "_zscore.csv");
		File zscoreCorrect = new File(myProt.getSource().getAbsolutePath().substring(0, myProt.getSource()
				.getAbsolutePath().indexOf(".pdb")) + "_zscoreCorrect.csv");

		File tempCSVfolder = ScoringGeneralHelpers.makeFolder(new File(myProt.getSource().getParent() + File
				.separator + "tempCSVs"));

		writeMatrixToCSV(resultCSV, intensityValueMatrix);
		writeMatrixToCSV(zscoreCSV, zValueMatrix);
		writeMatrixToCSV(zscoreCorrect, originalAcidsScore);

		writeTrueValueCSVs(zscoreCorrect, originalAcidsScore, originalPos);


	}

	private void writeTrueValueCSVs(File outputCSV, double[] originalAcidsScore, Integer[] originalPos) throws
			IOException {
		FileWriter FW = new FileWriter(outputCSV);

		String row = "";
		for (int j = 0; j < originalAcidsScore.length; j++) {
			row += originalPos[j] + ", " + originalAcidsScore[j] + "\n";
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

	private Integer[] getOriginalPositions(SimpleProtein myProt) throws InvalidPropertiesFormatException {
		List<Integer> positionArray = new ArrayList<>();
		for (SimpleProtein.ProtChain chain : myProt) {
			for (AminoAcid acid : chain) {
				positionArray.add(acidToIndex(acid.getName()));
			}
		}
		return positionArray.toArray(new Integer[positionArray.size()]);
	}


}
