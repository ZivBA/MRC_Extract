package utils.Scoring;

import utils.fileUtilities.MRC_Map_New;
import utils.molecularElements.AminoAcid;
import utils.molecularElements.SimpleAtom;
import utils.molecularElements.SimpleProtein;

import java.io.File;
import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.InvalidPropertiesFormatException;
import java.util.List;

/**
 * Created by zivben on 09/08/15.
 */
public class MRC_Score {
	private MRC_Map_New myMap;
	private SimpleProtein myProt;
	private double[][] resultMatrix;
	private double[][] zValueGrid;

	public MRC_Score(MRC_Map_New myMap, SimpleProtein myProt) {
		this.myMap = myMap;
		this.myProt = myProt;
		resultMatrix = new double[20][myProt.getLegnth()];

	}

	public MRC_Score(String mapPath, String protPath) throws IOException {
		this(new MRC_Map_New(mapPath), new SimpleProtein(new File(protPath)));
	}

	public double[][] scoreProtein() throws IOException {
		myProt.createPermutations();

		processSCWRLfolder(myProt.getProcessingFolder());

		return resultMatrix;
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

		for (SimpleProtein.ProtChain chain : tempProt) {
			for (AminoAcid residue : chain) {
				double scoreSum = 0;
				for (SimpleAtom atom : residue) {
					scoreSum += scoreSingleAtom(atom);
				}
				// put the residue cumulative score in the respective position in the resultMatrix
				// 1st position is the residue number minus the number of first residue, normalizing the
				// sequence number to a zero-start position. 2nd position is the index by acid name.
				resultMatrix[acidToIndex(residue.getName())][residue.getSeqNum() - tempProt.getSequenceBias()
						] = scoreSum;
			}
		}
	}

	private double scoreSingleAtom(SimpleAtom atom) {
		float[] coords = atom.getAtomCoords();

		atom.setAtomScore(myMap.val(coords[0], coords[1], coords[2]));
		return atom.getAtomScore();
	}
	
	
	public void calcZvalue() {
		double tempAvg[] = new double[20];
		double tempStD[] = new double[20];
		zValueGrid = new double[resultMatrix.length][resultMatrix[0].length];
		int counter;
		for (int i = 0; i < resultMatrix.length; i++) {
			for (int j = 0; j < resultMatrix[i].length; j++) {
				tempAvg[i] += resultMatrix[i][j];
			}
			tempAvg[i] = tempAvg[i] / resultMatrix[i].length;
			for (int j = 0; j < resultMatrix[i].length; j++) {
				tempStD[i] += Math.pow(resultMatrix[i][j] - tempAvg[i], 2);
			}
			tempStD[i] = Math.sqrt(tempStD[i] / resultMatrix[i].length);
		}
		
		for (int i = 0; i < resultMatrix.length; i++) {
			for (int j = 0; j < resultMatrix[i].length; j++) {
				zValueGrid[i][j] = (resultMatrix[i][j] - tempAvg[i]) / tempStD[i];
			}
		}
		
		
	}

	public void dispHist() {


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
