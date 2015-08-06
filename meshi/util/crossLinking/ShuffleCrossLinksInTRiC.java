package meshi.util.crossLinking;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.TRIC.PutUnitInAnyTopPosition;
import meshi.util.file.File2StringArray;


/**
 * In TRiC I enumerated the different positions as follows (16 in all):
 * {pos1 , pos2 , ... , pos8 , pos-1 , pos-2 , ... , pos-8}. 
 *
 */
public class ShuffleCrossLinksInTRiC {
	
	protected double[][][] activePositions = null;
	// These are the STRUCTURED active sites that are part of an INTER-unit xl. 
	protected String[] sitesStrings = {"A",
			"A",
			"A",
			"A",
			"A",
			"A",
			"B",
			"B",
			"B",
			"B",
			"B",
			"B",
			"B",
			"B",
			"D",
			"D",
			"D",
			"D",
			"D",
			"D",
			"E",
			"E",
			"E",
			"E",
			"G",
			"G",
			"G",
			"G",
			"H",
			"H",
			"Q",
			"Q",
			"Q",
			"Q",
			"Q",
			"Q",
			"Z",
			"Z"	};
	protected int[] sitesResNums = {33,
			109,
			126,
			272,
			494,
			499,
			119,
			120,
			135,
			176,
			230,
			284,
			431,
			473,
			32,
			58,
			68,
			129,
			387,
			398,
			42,
			49,
			279,
			496,
			21,
			78,
			249,
			427,
			55,
			109,
			37,
			326,
			400,
			439,
			459,
			466,
			199,
			222};

	public ShuffleCrossLinksInTRiC() {
		activePositions = new double[sitesResNums.length][16][3];
	}

	public int[] shuffle(int numberOfxl , double violCutoff, int numberOfShuffles) {
		String units = "ABGDEHQZ";
		String[] ringPerms = File2StringArray.f2a("ringPermute.txt");
		int[] histogram = new int[numberOfxl+1];
		
		for (int shuf = 0 ; shuf<numberOfShuffles ; shuf++) {
			if (shuf%1000==0) {
				System.out.print("Doing: " + shuf);
			}
			// Make the shuffled XLs
			int[] xlAS1 = new int[numberOfxl];
			int[] xlAS2 = new int[numberOfxl];
			for (int c=0 ; c<numberOfxl ; c++) {
				xlAS1[c] = (int) Math.floor(Math.random()*sitesResNums.length);
				do {
					xlAS2[c] = (int) Math.floor(Math.random()*sitesResNums.length);				
				} while (sitesStrings[xlAS2[c]].equals(sitesStrings[xlAS1[c]]));
			}
			//			// injecting the true XLs
			//		CrosslinkVector xlVec = (new CrosslinkVector("TRIC_ATP_xlinks_result.txt")).filterOutIntraUnit();
			//		String upperSeqq = "AGZQHEBD";
			//		String lowerSeqq = "HEBDAGZQ";
			//		AtomList fullComplex = PutUnitInAnyTopPosition.buildFullComplex(upperSeqq, lowerSeqq);
			//		int counter = 0;
			//		for (Crosslink xl : xlVec) {
			//            if ((fullComplex.findAtomInList("CA", xl.protName1(), xl.absPos1())!=-1) &&
			//            	(fullComplex.findAtomInList("CA", xl.protName2(), xl.absPos2())!=-1)) {
			//            	int indexHere = find(xl.protName1(), xl.absPos1());
			//            	xlAS1[counter] = indexHere;
			//            	indexHere = find(xl.protName2(), xl.absPos2());
			//            	xlAS2[counter] = indexHere;
			//    			counter++;
			//            }
			//        }	

			// Running the 8! permutations
			for (int perm = 0 ; perm<ringPerms.length ; perm++) {
				for (int reg = 0 ; reg<8 ; reg++) {
					String upperSeq = "";
					String lowerSeq = "";
					for (int c=0 ; c<8 ; c++) {
						int permAtCup = Integer.parseInt(""+ringPerms[perm].charAt(c))-1;
						int permAtCdown = Integer.parseInt(""+ringPerms[perm].charAt((reg + c) % 8))-1;
						upperSeq += units.charAt(permAtCup);
						lowerSeq += units.charAt(permAtCdown);
					}

					// Do analysis here on one of the 8! pemutaions.
					int consistCounter = 0;
					for (int xlC=0 ; xlC<xlAS1.length ; xlC++) {
						int indInUpperStringOf1 = upperSeq.indexOf(sitesStrings[xlAS1[xlC]].charAt(0));
						int indInLowerStringOf1 = lowerSeq.indexOf(sitesStrings[xlAS1[xlC]].charAt(0));
						int indInUpperStringOf2 = upperSeq.indexOf(sitesStrings[xlAS2[xlC]].charAt(0));
						int indInLowerStringOf2 = lowerSeq.indexOf(sitesStrings[xlAS2[xlC]].charAt(0));
						double minimalDistance = Double.MAX_VALUE;
						double dis = Double.MAX_VALUE;
						// in ring
						dis = norm(activePositions[xlAS1[xlC]][indInUpperStringOf1] , 
								activePositions[xlAS2[xlC]][indInUpperStringOf2]);
						if (dis<minimalDistance) {
							minimalDistance = dis;
						}
						// upper-lower
						dis = norm(activePositions[xlAS1[xlC]][indInUpperStringOf1] , 
								activePositions[xlAS2[xlC]][8+indInLowerStringOf2]);
						if (dis<minimalDistance) {
							minimalDistance = dis;
						}
						// lower-upper
						dis = norm(activePositions[xlAS1[xlC]][8+indInLowerStringOf1] , 
								activePositions[xlAS2[xlC]][indInUpperStringOf2]);
						if (dis<minimalDistance) {
							minimalDistance = dis;
						}
						if (minimalDistance<violCutoff) {
							consistCounter++;
						}
					}
					histogram[consistCounter]++;
				}	// of registration loop
			}	// of permutaion loop	
		}	// of shuffle loop	
		
		return histogram;		
	}
	
	/* 
	 * Returns the index in 'sitesResNums' where the desired active site is found. 
	 * Otherwise, returns -1 if the desired active site is not found.
	 */
	public int find(String chainID, int resNum) {
		for (int c=0 ; c<sitesResNums.length ; c++) {
			if ((sitesResNums[c]==resNum) && (sitesStrings[c].equals(chainID))) {
				return c;
			}
		}
		return -1;
	}
	
	/*
	 * This TRiC specific coding will register the (x,y,z) CB of a certain residue in one of the 16 positions.
	 */
	public void registerOneActiveSiteInOneComplex(String chainID, int resNum, int complexPos, int activeSiteIndex) {
		int posInM8to8 = 0;
		if (complexPos<8) {
			posInM8to8 = complexPos+1;
		}
		else {
			posInM8to8 = 7 - complexPos;
		}
		AtomList inPos = PutUnitInAnyTopPosition.putAnyUnitInAnyPosition(chainID, posInM8to8);
		Atom atom = inPos.findAtomInList("CB", resNum);
		activePositions[activeSiteIndex][complexPos][0] = atom.x();
		activePositions[activeSiteIndex][complexPos][1] = atom.y();
		activePositions[activeSiteIndex][complexPos][2] = atom.z();
	}
	
	/*
	 * This TRiC specific coding will register an active site in all 16 positions.
	 */
	public void registerOneActiveSiteInAllPositions(int activeSiteIndex) {
		String activeSiteString = sitesStrings[activeSiteIndex];
		int activeSiteResNum = sitesResNums[activeSiteIndex];
		for (int c=0 ; c<16 ; c++) {
			registerOneActiveSiteInOneComplex(activeSiteString, activeSiteResNum, c, activeSiteIndex);
		}
	}
	
	/*
	 * This TRiC specific coding will register an all active sites in all 16 positions.
	 */
	public void registerAllActiveSites() {
		for (int c=0 ; c<sitesStrings.length ; c++) {
			registerOneActiveSiteInAllPositions(c);			
		}
	}
	
	/**
	 * Norm of 2 3-tupples
	 */
	private static double norm(double[] v1, double[] v2) {
		return Math.sqrt((v1[0]-v2[0])*(v1[0]-v2[0]) +
				(v1[1]-v2[1])*(v1[1]-v2[1]) +
				(v1[2]-v2[2])*(v1[2]-v2[2]));
	}
		

	public static void main(String[] args) {
		ShuffleCrossLinksInTRiC shuf = new ShuffleCrossLinksInTRiC();
		shuf.registerAllActiveSites();
		System.out.println("Starting shuffling NOW!!");
		int [] hist = shuf.shuffle(18,28.00,100000);
		for (int c=0 ; c<hist.length ; c++) {
			System.out.println(c + "  " + hist[c]);
		}
	}
	
}
