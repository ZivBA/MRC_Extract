package meshi.applications.connexin;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.StringTokenizer;

import meshi.util.file.File2StringArray;

public class CreateFullPDBs {
	
	public static void main(String[] args) throws Exception {
		String fileName= args[0]; // The data file
		String outputFileNamePrefix = args[1]; // File name prefix
		String[] dataFile = File2StringArray.f2a(fileName);
		String[] masterFile = File2StringArray.f2a(dataFile[0]);
		int lineInDataFile = 1;
		boolean endOfFile = false;
		while (!endOfFile) {
			String resType = dataFile[lineInDataFile++];
			String resNum = dataFile[lineInDataFile++];
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outputFileNamePrefix + "_" + resType + "_" + resNum + ".pdb")));
			int atomCounterToWritten=1;
			int lineInMasterFile = 0;
			for (String nextLine = dataFile[lineInDataFile++] ; !nextLine.equals("TER") ; ) {
				int atomNumberInDataLine = Integer.parseInt(nextLine.substring(6,11).trim());
				for (boolean atomNumbersMatched=false ; (lineInMasterFile<masterFile.length) & (!atomNumbersMatched) ; lineInMasterFile++) {
					if (masterFile[lineInMasterFile].startsWith("ATOM") | masterFile[lineInMasterFile].startsWith("HETATM")) {
						if (atomCounterToWritten == atomNumberInDataLine) {
							atomNumbersMatched = true;
							lineInMasterFile += 4;
							String chain = nextLine.substring(21,22);
							for (boolean takeLine = true ; takeLine ; ) {
								pw.println(nextLine.substring(0,6)+printWithWith5Spaces(atomCounterToWritten)+nextLine.substring(11));
								//System.out.println(nextLine.substring(0,6)+printWithWith5Spaces(atomCounterToWritten)+nextLine.substring(11));
								if (!nextLine.equals("TER") && !nextLine.substring(0,6).equals("REMARK")) {
									atomCounterToWritten++;
								}
								nextLine = dataFile[lineInDataFile++];
								if (nextLine.equals("TER") || !nextLine.substring(21,22).equals(chain)) {
									takeLine = false;
								}
							}
						}
						else {
							pw.println(masterFile[lineInMasterFile].substring(0,6)+printWithWith5Spaces(atomCounterToWritten)+masterFile[lineInMasterFile].substring(11));
							//System.out.println(masterFile[lineInMasterFile].substring(0,6)+printWithWith5Spaces(atomCounterToWritten)+masterFile[lineInMasterFile].substring(11));							
							atomCounterToWritten++;
						}
					}
					else {
						if (masterFile[lineInMasterFile].startsWith("TER")) {
							pw.println("TER   "+printWithWith5Spaces(atomCounterToWritten)+"   ");
							//System.out.println("TER   "+printWithWith5Spaces(atomCounterToWritten)+"   ");
							atomCounterToWritten++;
						}
						else {
							pw.println(masterFile[lineInMasterFile]);
							//System.out.println(masterFile[atomCounterToWritten]);
						}
					}
				}
			}
			// "TER" occurred in data file.
			for ( ; lineInMasterFile<masterFile.length ; lineInMasterFile++) {
				if (masterFile[lineInMasterFile].startsWith("ATOM") | masterFile[lineInMasterFile].startsWith("HETATM")) {
					pw.println(masterFile[lineInMasterFile].substring(0,6)+printWithWith5Spaces(atomCounterToWritten)+masterFile[lineInMasterFile].substring(11));
					//System.out.println(masterFile[lineInMasterFile].substring(0,6)+printWithWith5Spaces(atomCounterToWritten)+masterFile[lineInMasterFile].substring(11));							
					atomCounterToWritten++;
				}
				else {
					if (masterFile[lineInMasterFile].startsWith("TER")) {
						pw.println("TER   "+printWithWith5Spaces(atomCounterToWritten)+"   ");
						//System.out.println("TER   "+printWithWith5Spaces(atomCounterToWritten)+"   ");
						atomCounterToWritten++;
					}
					else {
						pw.println(masterFile[lineInMasterFile]);
						//System.out.println(masterFile[lineInMasterFile]);
					}
				}
			}
			if (lineInDataFile==dataFile.length) {
				endOfFile = true;
			}
			pw.close();
		}
	}
	
	static String printWithWith5Spaces(int num) {
		if (num>9999)
			return ""+num;
		else if (num>999)
			return " "+num;
		else if (num>99)
			return "  "+num;
		else if (num>9)
			return "   "+num;
		else 
			return "    "+num;
	}
	
}
