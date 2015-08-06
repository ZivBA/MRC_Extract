package meshi.util.crossLinking;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.Vector;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.dssp.DSSP;
import meshi.util.file.File2StringArray;
import meshi.util.file.MeshiWriter;

public abstract class CrosslinkVector extends Vector<Crosslink> {
	
	private static final long serialVersionUID = 1L;
	private String creatingFile = null;
	private double violationCutoff = 27.55; //Angs
	
	/**
	 * An empty Crosslining vector, with 'null' 'creatingFile'.
	 */
	public CrosslinkVector() {}
	
	/**
	 * Building an vector of XL objects from the Abersold lab excel file (fileType=0), or from the Kalisman lab excel file (fileType=1).
	 */
	public CrosslinkVector(String filename , int fileType) {
		this.creatingFile = filename;
		String[] lines = File2StringArray.f2a(filename);
		for (String lineStr : lines) {
			add(new Crosslink(lineStr,fileType));
			lastElement().set_MS_file_name(filename);
		}		
	}

	
	public String toString() {
		String outString = "";
		for (int c=0; c<size() ; c++) {
			outString += ((c+1) + ".   " + get(c).toString());
		}
		return outString;
	}
	
	
	public String creatingFile() {
		return creatingFile;
	}
	
	public void setCreatingFile(String str) {
		creatingFile = str;
	}	
	
	public double violationCutoff() {
		return violationCutoff;
	}
	
	public void setViolationCutoff(double newVal) {
		violationCutoff = newVal;
	}
	
	/**
	 * Will appand the vector 'xlVec' at the end of 'this'.
	 */
	public void addVec(CrosslinkVector xlVec) {
		for (Crosslink xl : xlVec) {
			add(xl);
		}
	}

	/**
	 * Will appand the all the enteries of the vector 'xlVec' 
	 * that are not already in 'this' vector at the end of 'this'.
	 */
	public void addVecNonRedundant(CrosslinkVector xlVec) {
		for (Crosslink xl : xlVec) {
			if (find(xl.protName1(), xl.absPos1(), xl.protName2(), xl.absPos2()) == null) {
				add(xl);
			}
		}
	}

	/**
	 * this method will filter out redundancies of equal positions (not peptides). 
	 */
	public CrosslinkVector filterOutRedundancies() {
		CrosslinkVector newVec = createEmptyVector();
		for (Crosslink xl : this) {
			if (newVec.find(xl.protName1(), xl.absPos1(), xl.protName2(), xl.absPos2()) == null) {
				newVec.add(xl);
			} 
		}
		return newVec;
	}
	
	
	/**
	 * Will return a list without cross-links within the same protein. 
	 * Note: that this could also remove inter-protein links between homo-dimers. 
	 */
	public CrosslinkVector filterOutIntraUnit() {
		CrosslinkVector newVec = createEmptyVector();
		for (Crosslink xl : this) {
			if (!xl.protName1().equals(xl.protName2())) {
				newVec.add(xl);
			} else if (xl.absPos1() == xl.absPos2()) {
				newVec.add(xl);
			}
		}
		return newVec;
	}
	
	/**
	 * Will return only cross-links within the same protein, or between homo-dimers. 
	 */
	public CrosslinkVector filterOutInterUnit() {
		CrosslinkVector newVec = createEmptyVector();
		for (Crosslink xl : this) {
			if (xl.protName1().equals(xl.protName2())) {
				newVec.add(xl);
			}
		}
		return newVec;
	}
		
	/**
	 * Will return only cross-links above the threshold score. 
	 */
	public CrosslinkVector filterAboveScore(double thresholdScore) {
		CrosslinkVector newVec = createEmptyVector();
		for (Crosslink xl : this) {
			if (xl.score()>=thresholdScore) {
				newVec.add(xl);
			}
		}
		return newVec;
	}
	
	/**
	 * Will return only cross-links that stem from structured regions in the model. 
	 */
	public CrosslinkVector filterStructured(AtomList fullComplex) {
		CrosslinkVector newVec = createEmptyVector();
		for (Crosslink xl : this) {
			if ((fullComplex.findAtomInList("CA", xl.protName1(), xl.absPos1())!=-1) && 
				(fullComplex.findAtomInList("CA", xl.protName2(), xl.absPos2())!=-1)) {
				newVec.add(xl);
			}
		}
		return newVec;
	}

	
	/**
	 * Will return only cross-links where at least one of the residues is of protein that starts with a letter in the string.
	 **/
	public CrosslinkVector filterProteinsInSet(String proteinNames) {
		CrosslinkVector newVec = createEmptyVector();
		for (Crosslink xl : this) {
			if ( (proteinNames.indexOf(xl.protName1().charAt(0)) != -1) | (proteinNames.indexOf(xl.protName2().charAt(0)) != -1)) {
				newVec.add(xl);
			}
		}
		return newVec;
	}
	

	/**
	 * Will remove cross-links where at least one of the residues is of protein that starts with a letter in the string.
	 **/
	public CrosslinkVector filterOutProteinsInSet(String proteinNames) {
		CrosslinkVector newVec = createEmptyVector();
		for (Crosslink xl : this) {
			if ( (proteinNames.indexOf(xl.protName1().charAt(0)) == -1) & (proteinNames.indexOf(xl.protName2().charAt(0)) == -1)) {
				newVec.add(xl);
			}
		}
		return newVec;
	}
	
	/**
	 * Will remove cross-links where both residues are of proteins that starts with a letter in the string.
	 **/
	public CrosslinkVector filterOutBothProteinsInSet(String proteinNames) {
		CrosslinkVector newVec = createEmptyVector();
		for (Crosslink xl : this) {
			if ( (proteinNames.indexOf(xl.protName1().charAt(0)) == -1) | (proteinNames.indexOf(xl.protName2().charAt(0)) == -1)) {
				newVec.add(xl);
			}
		}
		return newVec;
	}

	
	
	public boolean isResidueCrossLinked(String protName, int resNum) {
		for (Crosslink xl : this) {
			if (xl.protName1().equals(protName) && 
					(xl.absPos1()==resNum)) {
				return true;
			}
			if (xl.protName1().equals(protName) && 
					(xl.absPos2()==resNum)) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * This method is TRiC specific.
	 * This method is Lysine-specific.
	 */
	public void printSASAstats(DSSP dssp, AtomList fullStruct) {
		for (int c=0 ; c<fullStruct.size() ; c++) {
			if (fullStruct.atomAt(c).name().equals("CB") && 
					fullStruct.atomAt(c).residueName().equals("LYS")) {
				Atom atom = fullStruct.atomAt(c);
				if (atom.chain().equals("A") || 
						atom.chain().equals("B") ||
						atom.chain().equals("G") ||
						atom.chain().equals("D") ||
						atom.chain().equals("E") ||
						atom.chain().equals("H") ||
						atom.chain().equals("Q") ||
						atom.chain().equals("Z")) {
					if (isResidueCrossLinked(atom.chain(), atom.residueNumber())) {
						System.out.println(/*atom.chain() + " " + atom.residueNumber() + " " + */
								dssp.ACCofRes(atom.residueNumber(), atom.chain().charAt(0)) + " " + 
								dssp.relACCofRes(atom.residueNumber(), atom.chain().charAt(0)) + " 1");
					}
					else {
						System.out.println(/*atom.chain() + " " + atom.residueNumber() + " " + */ 
								dssp.ACCofRes(atom.residueNumber(), atom.chain().charAt(0)) + " " + 
								dssp.relACCofRes(atom.residueNumber(), atom.chain().charAt(0)) + " 0");						
					}					
				}				
			}
		}
	}
	
	/**
	 * This will test the consistency/inconsist of the XL vector with the fullStruct.
	 * If checkProt1 is "X" then all the vector is checked. If not, then only cross-linking between 
	 * the proteins in 'checkProt1' and 'checkProt2' are checked. The results are in the the form of 
	 * [consistencies , unstructured , violations]
	 * This method is Lysine-specific.
	 */
	public int[] checkConsistencyWithStructure(AtomList fullStruct, String checkProt1, String checkProt2 , boolean toPrint) {
		int[] result = {0,0,0};
		AtomList struct = fullStruct.filter(new AtomList.KCA_Filter());
		int consistentCounter = 0;
		int terminiCounter = 0;
		int inConsistentCounter = 0;
		for (int c=0 ; c<size() ; c++) {
			Crosslink xl = get(c);
			if (checkProt1.equals("X") || 
					(xl.protName1().equals(checkProt1) && xl.protName2().equals(checkProt2)) ||
					(xl.protName1().equals(checkProt2) && xl.protName2().equals(checkProt1))) {
				String chains1 = identicalChains(xl.protName1()); 
				String chains2 = identicalChains(xl.protName2());
				if (xl.protName1().equals(xl.protName2()) && (xl.absPos1()==xl.absPos2())) { // This test will also check sure XLs across homo-dimers.
					chains1 = chains1.substring(0,1);
					chains2 = chains2.substring(1);
				}
				int resNum1 = xl.absPos1();
				int resNum2 = xl.absPos2();
				boolean consist = false;
				boolean termini = false;
				double minimalDistance = Double.MAX_VALUE;
				for (int chainInd1=0 ; chainInd1<chains1.length() ; chainInd1++) {
					for (int chainInd2=0 ; chainInd2<chains2.length() ; chainInd2++) {
						Atom a1 = struct.findAtomInListReturningAtom("CA", ""+chains1.charAt(chainInd1), resNum1);
						Atom a2 = struct.findAtomInListReturningAtom("CA", ""+chains2.charAt(chainInd2), resNum2);
						if ((a1!=null) && (a2!=null) && twoAtomsConsistent(a1,a2)) {
							consist = true;
							if (a1.distanceFrom(a2)<minimalDistance) {
								minimalDistance = a1.distanceFrom(a2);
							}								
						}
						else if ((a1==null) || (a2==null)) {
							termini = true;
						}
						else {
							if (a1.distanceFrom(a2)<minimalDistance) {
								minimalDistance = a1.distanceFrom(a2);
							}
						}
					}
				}
				if (minimalDistance < Double.MAX_VALUE) {
//				double rounded = Math.round(minimalDistance);
//				System.out.print(" " + rounded);
				}
				if (consist) {
					consistentCounter++;
					result[0]++;
				}
				else if (termini) {
					terminiCounter++;
					result[1]++;
					if (toPrint) {
						System.out.println("Termini:\n" + xl);
					}
				}
				else {
					result[2]++;
					inConsistentCounter++;
					if (toPrint) {
						System.out.println("Inconsistent:\nD" + minimalDistance + "D\n" + xl);
					}
				}
			}
		}
		if (toPrint) {
			System.out.println("\n-------------------------------------------\n" + 
				"Total:    good=" + consistentCounter + "    termini=" + terminiCounter + "    bad=" + inConsistentCounter);
		}
		return result;
	}


	/**
	 * This will return a string with the distances of the cross-links on the model. If a cross-linked residue is missing the nearest residue
	 * up to 'missingResTol' will be assigned. If any of the cross-linked residues are not in the model, NA will be printed. 
	 */
	public String printDistancesOnStructure(AtomList modelFull, int missingResTol) {
		AtomList model = modelFull.CAFilter();
		String outString = "";
		for (Crosslink xl : this) {
			// Looking for closest to residues
			Atom nearest1 = null;
			int resDis1 = Integer.MAX_VALUE;
			Atom nearest2 = null;
			int resDis2 = Integer.MAX_VALUE;
			for (int c=0 ; c<model.size() ; c++) {
				if (model.atomAt(c).chain().equals(xl.protName1())) {
					if (Math.abs(xl.absPos1()-model.atomAt(c).residueNumber())<resDis1) {
						resDis1 = Math.abs(xl.absPos1()-model.atomAt(c).residueNumber());
						nearest1 = model.atomAt(c);
					}
				}
				if (model.atomAt(c).chain().equals(xl.protName2())) {
					if (Math.abs(xl.absPos2()-model.atomAt(c).residueNumber())<resDis2) {
						resDis2 = Math.abs(xl.absPos2()-model.atomAt(c).residueNumber());
						nearest2 = model.atomAt(c);
					}
				}
			}
			
//			if ((nearest1==null) || (nearest2==null) || (resDis1>missingResTol) || (resDis2>missingResTol)) {
//				outString += "NA\n";
//			}
//			else {
//				outString += ((((int) (10*nearest1.distanceFrom(nearest2)))/10.0) + "\n");
//			}

			if ((nearest1==null) || (nearest2==null)) {
				outString += "NA  NA\n";
			}
			else if ((resDis1<missingResTol) & (nearest1.occupancy()<2) & (resDis2<missingResTol) & (nearest2.occupancy()<2)) {
				if (nearest1.distanceFrom(nearest2)>32.0) 
					outString += ((((int) (10*nearest1.distanceFrom(nearest2)))/10.0) + "^^^  NA\n");
				else
					outString += ((((int) (10*nearest1.distanceFrom(nearest2)))/10.0) + "  NA\n");
			}
			else if ((nearest1.occupancy()>1) | (nearest2.occupancy()>1)) {
				if (nearest1.distanceFrom(nearest2)==0.0)
					outString += ("NA  NA\n");
				else
					outString += ("NA  " + (((int) (10*nearest1.distanceFrom(nearest2)))/10.0) + "\n");
			} 
			else
				outString += "NA  NA\n";;
		
		}		
		return outString;
	}
	
	
	
	/**
	 * This method is similar to 'checkConsistencyWithStructure', but gives the sum of violations score.
	 * If checkProt1 is "X" then all the vector is checked. If not, then only cross-linking between 
	 * the proteins in 'checkProt1' and 'checkProt2' are checked. The results are in the the form of 
	 * [violations-24.0Angs, violations-28.0Angs , sumViol-24.0Angs , sumViol-28.0Angs]
	 * This method is Lysine-specific.
	 */
	public double[] checkConsistencyWithStructureSumOfViolations(AtomList fullStruct, String checkProt1, String checkProt2 , boolean toPrint) {
		double[] result = {0.0 , 0.0 , 0.0 , 0.0};
		AtomList struct = fullStruct.filter(new AtomList.KCA_Filter()); // Cange here to CB
		for (int c=0 ; c<size() ; c++) {
			Crosslink xl = get(c);
			if (checkProt1.equals("X") || 
					(xl.protName1().equals(checkProt1) && xl.protName2().equals(checkProt2)) ||
					(xl.protName1().equals(checkProt2) && xl.protName2().equals(checkProt1))) {
				String chains1 = identicalChains(xl.protName1()); 
				String chains2 = identicalChains(xl.protName2());
				if (xl.protName1().equals(xl.protName2()) && (xl.absPos1()==xl.absPos2())) { // This test will also check sure XLs across homo-dimers.
					chains1 = chains1.substring(0,1);
					chains2 = chains2.substring(1);
				}
				int resNum1 = xl.absPos1();
				int resNum2 = xl.absPos2();
				boolean consist = false;
				boolean termini = false;
				double minimalDistance = Double.MAX_VALUE;
				for (int chainInd1=0 ; chainInd1<chains1.length() ; chainInd1++) {
					for (int chainInd2=0 ; chainInd2<chains2.length() ; chainInd2++) {
						Atom a1 = struct.findAtomInListReturningAtom("CA", ""+chains1.charAt(chainInd1), resNum1);  // Cange here to CB
						Atom a2 = struct.findAtomInListReturningAtom("CA", ""+chains2.charAt(chainInd2), resNum2);  // Cange here to CB
						if ((a1!=null) && (a2!=null) && twoAtomsConsistent(a1,a2)) {
							consist = true;
							if (a1.distanceFrom(a2)<minimalDistance) {
								minimalDistance = a1.distanceFrom(a2);
							}								
						}
						else if ((a1==null) || (a2==null)) {
							termini = true;
						}
						else {
							if (a1.distanceFrom(a2)<minimalDistance) {
								minimalDistance = a1.distanceFrom(a2);
							}
						}
					}
				}
				if (minimalDistance < Double.MAX_VALUE) {
//				double rounded = Math.round(minimalDistance);
//				System.out.print(" " + rounded);
				}
				if (consist) {
//
				}
				else if (termini) {
					if (toPrint) {
						System.out.println("Termini:\n" + xl);
					}
				}
				else {
					result[0]++;
					result[2] += (minimalDistance-24.0);
					if (minimalDistance>28.0) {
						result[1]++;
						result[3] += (minimalDistance-28.0);					
					}
					if (toPrint) {
						System.out.println("Inconsistent:\nD" + minimalDistance + "D\n" + xl);
					}
				}
			}
		}
		if (toPrint) {
			System.out.println("\n-------------------------------------------\n" + 
				"Total:   #Viol24=" + result[0] + "    #Viol28=" + result[1] + "    SumViol24=" + result[2] + "    sumViol28=" + result[3]);
		}
		return result;
	}
	
	
	public boolean twoAtomsConsistent(Atom a1, Atom a2) {
		if (a1.distanceFrom(a2) < violationCutoff) {
			return true;
		}
		else {
			return false;
		}
	}
	
	
	abstract String identicalChains(String protName);

	abstract CrosslinkVector createEmptyVector();

	/**
	 * Check the XL data for errors when compared to the sequences in 'sequencesFile'. 
	 **/
	public void checkDataValidity(String sequencesFile) {
		boolean goodData = true;
		MySequenceList seqList = new MySequenceList(sequencesFile);
		for (Crosslink xl : this) {
			int[][] matches = seqList.findInSequences(xl.seq1());
			if (matches.length == 0) {
				System.out.println("Error: could not find a match in the sequences for this MS sequence:");
				System.out.println(xl.seq1() + "    in protein: " + xl.protName1() + "\n");
				goodData = false;
			} 
			else if (matches.length > 1) {
				System.out.println(matches.length + " Error: found more than 1 match for this MS sequence:");
				System.out.print(xl.seq1() + "    in protein: " + xl.protName1() + "        " + xl);
				System.out.println("These matches are in proteins (first 2 only):");
				System.out.println(seqList.get(matches[0][0]).title() + " " + matches[0][1]);
				System.out.println(seqList.get(matches[1][0]).title() + " " + matches[1][1] + "\n");
				goodData = false;
			} 
			else if ((xl.absPos1() < matches[0][1]) || (xl.absPos1() > (matches[0][1] + xl.seq1().length() -1))) {
				System.out.println("Error: The absolute position in this MS sequence is outside the sequence:");
				System.out.println(xl.seq1() + "    in protein: " + xl.protName1() + "     absolute position: " + xl.absPos1() + "\n");
				goodData = false;				
			}
			else if (seqList.get(matches[0][0]).seq().charAt(xl.absPos1()-1) != 'K') {
				System.out.println("Error: The absolute position in this MS sequence is not on a lysine:");
				System.out.println(xl.seq1() + "    in protein: " + xl.protName1() + "     absolute position: " + xl.absPos1() + "\n");
				goodData = false;								
			}
			matches = seqList.findInSequences(xl.seq2());
			if (matches.length == 0) {
				System.out.println("Error: could not find a match in the sequences for this MS sequence:");
				System.out.println(xl.seq2() + "    in protein: " + xl.protName2() + "\n");
				goodData = false;
			} 
			else if (matches.length > 1) {
				System.out.println(matches.length + " Error: found more than 1 match for this MS sequence:");
				System.out.print(xl.seq2() + "    in protein: " + xl.protName2() + "        " + xl);
				System.out.println("These matches are in proteins (first 2 only):");
				System.out.println(seqList.get(matches[0][0]).title() + " " + matches[0][1]);
				System.out.println(seqList.get(matches[1][0]).title() + " " + matches[1][1] + "\n");
				goodData = false;
			} 
			else if ((xl.absPos2() < matches[0][1]) || (xl.absPos2() > (matches[0][1] + xl.seq2().length() -1))) {
				System.out.println("Error: The absolute position in this MS sequence is outside the sequence:");
				System.out.println(xl.seq2() + "    in protein: " + xl.protName2() + "     absolute position: " + xl.absPos1() + "\n");
				goodData = false;				
			}
			else if (seqList.get(matches[0][0]).seq().charAt(xl.absPos2()-1) != 'K') {
				System.out.println("Error: The absolute position in this MS sequence is not on a lysine:");
				System.out.println(xl.seq2() + "    in protein: " + xl.protName2() + "     absolute position: " + xl.absPos1() + "\n");
				goodData = false;								
			}						
		}
		if (goodData) {
			System.out.println("Data is OK.");
		}
	}
		
	
	public String compareToOtherXLvec(CrosslinkVector otherVec) {
		String result = "These XLs are in: " + this.creatingFile() + " , but not in: " +
		otherVec.creatingFile() + "\n----------------------------------------------------------------------------\n";		
		// Checking one direction ->
		int counter = 0;
		for (Crosslink xl : this) {
			if (otherVec.find(xl.protName1(), xl.absPos1(), xl.protName2(), xl.absPos2())==null) {
				counter++;
				result += counter + ".  " +xl.toString(); 
			}
		}
		result += "\nThese XLs are in: " + otherVec.creatingFile() + " , but not in: " +
		this.creatingFile() + "\n----------------------------------------------------------------------------\n";
		// Checking the other direction <-
		counter = 0;
		for (Crosslink xl : otherVec) {
			if (this.find(xl.protName1(), xl.absPos1(), xl.protName2(), xl.absPos2())==null) {
				counter++;
				result += counter + ".  " +xl.toString(); 
				System.out.println(xl.absPos1()+"\n"+xl.absPos2());
			}
		}
		result += "\nThese XLs are shared:\n" + 
		"----------------------------------------------------------------------------\n";
		// Checking the other direction <-
		counter = 0;
		for (Crosslink xl : otherVec) {
			if (this.find(xl.protName1(), xl.absPos1(), xl.protName2(), xl.absPos2())!=null) {
				counter++;
				result += counter + ".  " +xl.toString();
			}
		}
		return result;
	}
	
    /**
     * Return a list of all the active (crosslinked) AND STRUCTURED residues in the vector. 
     */
    public String getActiveStructuredResidues(AtomList fullComlex) {
            String out = "";
            for (Crosslink xl : this) {
                    if (fullComlex.findAtomInList("CA", xl.protName1(), xl.absPos1())!=-1) {
                            out += xl.protName1() + " " + xl.absPos1() + "\n";
                    }
                    if (fullComlex.findAtomInList("CA", xl.protName2(), xl.absPos2())!=-1) {
                            out += xl.protName2() + " " + xl.absPos2() + "\n";
                    }
//                    out += xl.protName1() + " " + xl.absPos1() + "\n" + xl.protName2() + " " + xl.absPos2() + "\n";
            }
            return out;
    }
	
	
	/**
	 * Returns the first crosslink from the Vec with the properties in the parameters. If none is found, 
	 * returns null.
	 */
	public Crosslink find(String prot1, int pos1, String prot2, int pos2) {
		for (Crosslink xl : this) {
			if ((xl.absPos1()==pos1) && (xl.absPos2()==pos2) && xl.protName1().equals(prot1) 
					&& xl.protName2().equals(prot2)) {
				return xl;
			}
			if ((xl.absPos1()==pos2) && (xl.absPos2()==pos1) && xl.protName1().equals(prot2) 
					&& xl.protName2().equals(prot1)) {
				return xl;
			}
		}
		return null;
	}
	
	public String printListWithStructuralInfo(AtomList fullComplex) {
        String out = "";
		for (Crosslink xl : this) {
			String prot1StructStatus = null;
			String prot2StructStatus = null;
			if (fullComplex.findAtomInList("CA", xl.protName1(), xl.absPos1())!=-1) {
				prot1StructStatus = "S";
			}
			else {
				prot1StructStatus = "U";				
			}
	        if (fullComplex.findAtomInList("CA", xl.protName2(), xl.absPos2())!=-1) {
				prot2StructStatus = "S";
			}
			else {
				prot2StructStatus = "U";				
			}
	        out += xl.protName1() + " " + xl.absPos1() + " " + xl.protName2() + " " + xl.absPos2() + " " + prot1StructStatus + " " + prot2StructStatus + "\n";
		}
        return out;
    }
	
}
