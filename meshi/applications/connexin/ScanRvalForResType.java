package meshi.applications.connexin;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Vector;

import meshi.util.file.File2StringArray;

public class ScanRvalForResType {

	public static void main(String[] args) throws Exception {
		String pdbID = args[0]; // The PDB to work on
		String chainID = args[1]; // The Chain ID
		String identicalChains = args[2]; // The identical Chains
		String aaType = args[3]; // The aa type to scan
		
		String[] pristineMasterFile = File2StringArray.f2a(pdbID + "_" + chainID + ".MASTER.pdb");
		String[] resNumFile = File2StringArray.f2a(pdbID + "_" + chainID + ".RESNUM.pdb");
		PrintWriter outputWriter = new PrintWriter(new BufferedWriter(new FileWriter(pdbID + "_" + chainID + "_" + aaType + ".RVALS.pdb")));
		for (int resNumInd=0 ; resNumInd<resNumFile.length ; resNumInd++) {
			int resNum = Integer.parseInt(resNumFile[resNumInd]);
			String[] masterFile = new String[pristineMasterFile.length];  // New copy of master ALA file
			for (int c=0 ; c<masterFile.length ; c++) {
				masterFile[c] = new String(pristineMasterFile[c]);
			}
			for (int identChainInd=0 ; identChainInd<identicalChains.length() ; identChainInd++) {
				char workChain = identicalChains.charAt(identChainInd);
				int lineNumber = getAtomReturnLineNum(masterFile, workChain, resNum, "CA");
				if (lineNumber != -1) {
					lineNumber = getAtomReturnLineNum(masterFile, workChain, resNum, "N");
					masterFile[lineNumber] = masterFile[lineNumber].substring(0,17) + aaType + masterFile[lineNumber].substring(20);
					lineNumber = getAtomReturnLineNum(masterFile, workChain, resNum, "CA");
					masterFile[lineNumber] = masterFile[lineNumber].substring(0,17) + aaType + masterFile[lineNumber].substring(20);
					lineNumber = getAtomReturnLineNum(masterFile, workChain, resNum, "C");
					masterFile[lineNumber] = masterFile[lineNumber].substring(0,17) + aaType + masterFile[lineNumber].substring(20);
					lineNumber = getAtomReturnLineNum(masterFile, workChain, resNum, "O");
					masterFile[lineNumber] = masterFile[lineNumber].substring(0,17) + aaType + masterFile[lineNumber].substring(20);
					lineNumber = getAtomReturnLineNum(masterFile, workChain, resNum, "CB");
					if (lineNumber==-1) {
						System.out.println(workChain + " " + resNum + " ");
					}
					masterFile[lineNumber] = masterFile[lineNumber].substring(0,17) + aaType + masterFile[lineNumber].substring(20);
					if (aaType.equals("GLY")) {
						masterFile[lineNumber] = "REMARK  " + masterFile[lineNumber].substring(8);						
					}
				}
			}
			// running SCWRL
			printPDB(masterFile, "ReadyForSCWRL.pdb");
	        try {
	            Process p = Runtime.getRuntime().exec("Scwrl4 -i ReadyForSCWRL.pdb -o AfterSCWRL.pdb");
	            BufferedReader in = new BufferedReader( new InputStreamReader(p.getInputStream()) );
	            String line = null;
	            while ((line = in.readLine()) != null) {
	                System.out.println(line);
	            }
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
			System.out.println("SCWRL is finished");
			// Merging SCWRL and original file.
			String[] SCWRLfile = File2StringArray.f2a("AfterSCWRL.pdb");
			String[] readyPDB = mergeSCWRL(masterFile,SCWRLfile,false,true);
			printPDB(readyPDB, "ReadyForSFCHECK.pdb");
			// Running SFCHECK
			try {
	            Process p = Runtime.getRuntime().exec("sfcheck -f " + pdbID + "-sf.cif -m ReadyForSFCHECK.pdb -mem 450 -na 90000");
	            BufferedReader in = new BufferedReader( new InputStreamReader(p.getInputStream()) );
	            String line = null;
	            while ((line = in.readLine()) != null) {
	                System.out.println(line);
	            }
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
			// Extracting the R-values and correlation
			String[] RvalFile = File2StringArray.f2a("sfcheck.log");
			outputWriter.print(resNum + " ");
			boolean foundRval = false;
			for (int RvalFileCounter=0 ; RvalFileCounter<RvalFile.length ; RvalFileCounter++ ) {
				if (RvalFile[RvalFileCounter].length() > 21) {
					if (RvalFile[RvalFileCounter].substring(0, 21).equals(" R-factor           :") ) {
						outputWriter.print(Double.parseDouble(RvalFile[RvalFileCounter].substring(21)) + " ");
						foundRval = true;
					}
					if (RvalFile[RvalFileCounter].substring(0, 21).equals(" Correlation factor :") ) {
						outputWriter.println(Double.parseDouble(RvalFile[RvalFileCounter].substring(21)) + " ");
					}
				}
			}
			if (!foundRval) {
				outputWriter.println("\nCOURRPTED:  " + pdbID + "_" + chainID + "_" + aaType);
				outputWriter.close();
				throw new RuntimeException("Error in the run. Stopping.");
			}
		}
		outputWriter.close();
	}
	
	
	
	public static String[] mergeSCWRL(String[] REFpdbFile, String[] scwrlFile, boolean takeCBfromSCWRL, boolean renumber) {
		Vector<String> outPDB = new Vector<String>();
		int resNUM = -1;
		String TF = "-----------";
		int indSCWRL = 0;
		int atomNumberCounter = 0;
		for (int startC=0 ; startC<REFpdbFile.length ; startC++) {
			if (!REFpdbFile[startC].substring(0, 3).equals("ATO")) {
				String newLine = "";
	            if ( REFpdbFile[startC].substring(0, 6).equals("HETATM") | REFpdbFile[startC].substring(0, 3).equals("TER") ) {
	            	atomNumberCounter++;
	            	if (renumber) {
	            		newLine = REFpdbFile[startC].substring(0,6) + printWithWith5Spaces(atomNumberCounter) + REFpdbFile[startC].substring(11);
	            	}
	            	else {
	            		newLine = REFpdbFile[startC];
	            	}
	            }
	            else {
	            	newLine = REFpdbFile[startC]; 
	            }
				outPDB.add(newLine);
			}
			else {
				String newLine = "";
	            if ( REFpdbFile[startC].substring(0, 4).equals("ATOM") | REFpdbFile[startC].substring(0, 6).equals("HETATM") | REFpdbFile[startC].substring(0, 3).equals("TER") ) {
	            	atomNumberCounter++;
	            	if (renumber) {
	            		newLine = REFpdbFile[startC].substring(0,6) + printWithWith5Spaces(atomNumberCounter) + REFpdbFile[startC].substring(11);
	            	}
	            	else {
	            		newLine = REFpdbFile[startC];
	            	}
	            }
	            else {
	            	newLine = REFpdbFile[startC]; 
	            }
				outPDB.add(newLine);
				if (REFpdbFile[startC].substring(12,16).trim().equals("CA")) {
					resNUM = Integer.parseInt( REFpdbFile[startC].substring(22,26).trim() );
					TF = REFpdbFile[startC].substring(60,66);
				}
				String atomName = REFpdbFile[startC].substring(12,16).trim();
				if ( (atomName.equals("O") & takeCBfromSCWRL) | atomName.equals("CB") 	) { //  Take sidechain atoms from SCWRL file after the BB oxygen or the CB are written (depending where the CB is coming from)
					int maxLineToTake  = Math.min(indSCWRL+40 , scwrlFile.length);
					for (int searchIND = indSCWRL ; searchIND<maxLineToTake ; searchIND++) {
						if (scwrlFile[searchIND].substring(0, 3).equals("ATO")) {
							if ( Integer.parseInt( scwrlFile[searchIND].substring(22,26).trim() ) == resNUM) {
								atomName = scwrlFile[searchIND].substring(12,16).trim();
								if (atomName.equals("CA")) {
									indSCWRL = searchIND;
								}
		                        if (!atomName.equals("N") & !atomName.equals("CA") & !atomName.equals("C") & !atomName.equals("O") & 
		                                (!atomName.equals("CB") | takeCBfromSCWRL) &
		                                (atomName.charAt(0) != 'H') & 
		                                (atomName.charAt(0) != '1') & 
		                                (atomName.charAt(0) != '2') &
		                                (atomName.charAt(0) != '3')) { // Just SC atoms that are not hydrogens
		                        	atomNumberCounter++;
		        	            	if (renumber) {
		        	            		newLine = scwrlFile[searchIND].substring(0,6) + printWithWith5Spaces(atomNumberCounter) + scwrlFile[searchIND].substring(11,60) + TF + scwrlFile[searchIND].substring(66);
		        	            	}
		        	            	else {
		        	            		newLine = scwrlFile[searchIND].substring(0,60) + TF + scwrlFile[searchIND].substring(66);
		        	            	}		                        	
		                        	outPDB.add( newLine   );
		                        }	            
							}
						}
					}
				}
			}
		}
		return  outPDB.toArray(new String[0]);
	}

	
	public static int getAtomReturnLineNum(String[] file, char chain_id, int resNum, String atomName) {
		for (int c=0 ; c<file.length ; c++) {
			if (file[c].substring(0,3).equals("ATO")) {
				if (file[c].charAt(21)==chain_id) { 
					if (Integer.parseInt(file[c].substring(22,26).trim())==resNum) {
						if (file[c].substring(12,16).trim().equals(atomName)) {
							return c;
						}
					}
				}
			}
		}
		return -1;
	}
	
	public static void printPDB(String[] pdb, String fileName) throws Exception {
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(fileName)));
		for (int c=0 ; c<pdb.length ; c++) {
			pw.println(pdb[c]);			
		}
		pw.close();
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
