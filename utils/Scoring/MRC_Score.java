package utils.Scoring;

import utils.fileUtilities.FileProcessor;
import utils.fileUtilities.MRC_Map_New;
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

/**
 * Created by zivben on 09/08/15.
 */
public class MRC_Score {
	private MRC_Map_New myMap;
	private SimpleProtein myProt;
	private double[][] intensityValueMatrix;
	private double[][] zValueMatrix;
	private double[][] zValueCorrect;
	private double[] zValuesForOriginalAcids;

	public MRC_Score(MRC_Map_New myMap, SimpleProtein myProt) {
		this.myMap = myMap;
		this.myProt = myProt;
		intensityValueMatrix = new double[20][myProt.getLegnth()];

	}

	public MRC_Score(String mapPath, String protPath) throws IOException {
		this(new MRC_Map_New(mapPath), new SimpleProtein(new File(protPath)));
	}

	public static MRC_Score StartFromScratch(FileProcessor FP, String mrcpath) {

		try {
			ProteinActions.stripAndAllALAToFile(FP.getSource(), FP.getDest());
			SimpleProtein processedProt = new SimpleProtein(FP.getDest());

			MRC_Score scoreMap = new MRC_Score(new MRC_Map_New(mrcpath), processedProt);
			scoreMap.scoreProtein();
			scoreMap.calcZvalue();
			scoreMap.dispHist();
			return scoreMap;
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		}

	}

	public double[][] scoreProtein() throws IOException {
		myProt.createPermutations();

		processSCWRLfolder(myProt.getProcessingFolder());

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
					double scoreSum = 0;
					for (SimpleAtom atom : residue) {
						scoreSum += scoreSingleAtom(atom);
					}
					// put the residue cumulative score in the respective position in the
					// intensityValueMatrix
					// 1st position is the residue number minus the number of first residue, normalizing the
					// sequence number to a zero-start position. 2nd position is the index by acid name.
					intensityValueMatrix[testResidueIndex][testResiduePosition - tempProt.getSequenceBias()]
							= scoreSum;
				}
			}

		}


	}
	
	private double scoreSingleAtom(SimpleAtom atom) {
		float[] coords = atom.getAtomCoords();

		atom.setAtomScore(myMap.val(coords[0], coords[1], coords[2]));
		return atom.getAtomScore();
	}

	public void calcZvalue() throws InvalidPropertiesFormatException {
		double tempAvg[] = new double[20];
		double tempStD[] = new double[20];
		zValueMatrix = new double[intensityValueMatrix.length][intensityValueMatrix[0].length];
		zValuesForOriginalAcids = new double[intensityValueMatrix.length];
		int counter;
		for (int i = 0; i < intensityValueMatrix.length; i++) {
			for (int j = 0; j < intensityValueMatrix[i].length; j++) {
				tempAvg[i] += intensityValueMatrix[i][j];
			}
			tempAvg[i] = tempAvg[i] / intensityValueMatrix[i].length;
			for (int j = 0; j < intensityValueMatrix[i].length; j++) {
				tempStD[i] += Math.pow(intensityValueMatrix[i][j] - tempAvg[i], 2);
			}
			tempStD[i] = Math.sqrt(tempStD[i] / intensityValueMatrix[i].length);
		}

		for (int i = 0; i < intensityValueMatrix.length; i++) {
			for (int j = 0; j < intensityValueMatrix[i].length; j++) {
				zValueMatrix[i][j] = (intensityValueMatrix[i][j] - tempAvg[i]) / tempStD[i];
			}
		}

		Integer[] suspectedCorrectPositions = getOriginalPositions(myProt);

		for (int i = 0; i < suspectedCorrectPositions.length; i++) {
			zValuesForOriginalAcids[i] = zValueMatrix[suspectedCorrectPositions[i]][i];
		}

	}

	public void dispHist() throws IOException {

		File resultCSV = new File(myProt.getSource().getAbsolutePath().substring(0, myProt.getSource()
				.getAbsolutePath().indexOf(".pdb")) + "_resultMatrix.csv");
		File zscoreCSV = new File(myProt.getSource().getAbsolutePath().substring(0, myProt.getSource()
				.getAbsolutePath().indexOf(".pdb")) + "_zscore.csv");


		writeMatrixToCSV(resultCSV, intensityValueMatrix);
		writeMatrixToCSV(zscoreCSV, zValueMatrix);


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

	private int acidToIndex(String name) throws InvalidPropertiesFormatException {

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
