package meshi.util.crossLinking;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.TRIC.PutUnitInAnyTopPositionYeast;
import meshi.util.file.File2StringArray;
import meshi.util.file.MeshiWriter;

public class CrosslinkVectorTRiC extends CrosslinkVector  {

	private static final long serialVersionUID = 1L;

	public CrosslinkVectorTRiC() {
		// TODO Auto-generated constructor stub
	}

	public CrosslinkVectorTRiC(String filename , int fileType) {
		super(filename ,  fileType);
		// TODO Auto-generated constructor stub
	}

	/**
	 * A TRiC specific method that translate a chain ID from the top ring into the ID that I 
	 * gave it in the bottom ring.
	 */
	public String identicalChains(String protName) {
		String outS = "" + protName.charAt(0) + 
		PutUnitInAnyTopPositionYeast.matchChainLetterOnBottomRing(protName.charAt(0));  
		return outS;
	}

	public CrosslinkVector createEmptyVector() {
		return new CrosslinkVectorTRiC();
	}	
	
	public static void main(String[] args) {

		// Calculating the distances of the cross-links on the new HUJI data set of CCT. 
		AtomList fullComplex = new AtomList("G:\\CCT_HUJI\\4AOL.pdb");
		CrosslinkVectorTRiC xlVec = new CrosslinkVectorTRiC("G:\\CCT_HUJI\\Report_2015_02_18_04_CCT-XL_old_above055.txt", 1);
		for (Crosslink xl : xlVec) {
			String str = "N/A N/A";
			if ((fullComplex.findAtomInListReturningAtom("CA", xl.protName1(), xl.absPos1()) != null) &
					(fullComplex.findAtomInListReturningAtom("CA", ""+xl.protName2().charAt(0), xl.absPos2()) != null)) {
				double sameRingDis = fullComplex.findAtomInListReturningAtom("CA", xl.protName1(), xl.absPos1()).distanceFrom(
						fullComplex.findAtomInListReturningAtom("CA", ""+xl.protName2().charAt(0), xl.absPos2()) );
				double acrossRingDis = fullComplex.findAtomInListReturningAtom("CA", xl.protName1(), xl.absPos1()).distanceFrom(
						fullComplex.findAtomInListReturningAtom("CA", ""+(xl.protName2().toLowerCase().charAt(0)), xl.absPos2()+1000) );
				if ((sameRingDis>acrossRingDis) | (sameRingDis<0.01)) {
					if (acrossRingDis>30.0) {
						str = acrossRingDis + " 2Rings_viol";
					}
					else {
						str = acrossRingDis + " 2Rings";
					}
				}
				else {
					if (sameRingDis>30.0) {
						str = sameRingDis + " VIOL";
					}
					else {
						str = sameRingDis + " --";
					}					
				}
			}
			System.out.println(str);
		}				
		
		
//		CrosslinkVector xlVec = new CrosslinkVector("TRIC_ATP_xlinks_result.txt");
//		CrosslinkVector other = new CrosslinkVector("TRIC_xlinks_result.txt");
//		System.out.println(xlVec.compareToOtherXLvec(other));

		
//		xlVec.checkDataValidity("TRIC_sequences.txt");

		
//		// Checking intra-unit.
//		// IMPORTANT!! IMPORTANT!! IMPORTANT!! 
//		// Put this line: 'return "" + protName.charAt(0);' at the very begining of 'protName2Chains' or the results will be mixed-up.
//		String units = "ABGDEHQZ";
//		for (int c=0 ; c<units.length() ; c++) {
//			String singleChain = "" + units.charAt(c);
//			System.out.println("\n-----------------------------------------------------");
//			for (int qqq=0; qqq<30 ; qqq++) {
//				System.out.print(singleChain);
//			}
//			System.out.println();
//			AtomList al = new AtomList("C:\\Users\\Nir\\TriC\\HM_models_of_units\\1A6D\\Atemp_HM_"+singleChain+".pdb");
//			al.setChain(singleChain);
//			xlVec.checkConsistencyWithStructure(al, singleChain, singleChain, true);
//			System.out.println();
//		}

//		// Checking intra-ring pairs
//		String units = "ABGDEHQZ";
//		AtomList pos1 = new AtomList("C:\\Users\\Nir\\TriC\\HM_models_of_units\\1A6D\\pos_1.pdb");
//		AtomList pos2 = new AtomList("C:\\Users\\Nir\\TriC\\HM_models_of_units\\1A6D\\pos_2.pdb");
//		for (int c=0 ; c<units.length() ; c++) {
//			for (int d=0 ; d<units.length() ; d++) {
//				System.out.println("\n----------- " + units.charAt(c) + units.charAt(d) + " -----------");
//				AtomList alLeft = new AtomList("C:\\Users\\Nir\\TriC\\HM_models_of_units\\1A6D\\Atemp_HM_"+units.charAt(c)+".pdb");
//				AtomList alRight = new AtomList("C:\\Users\\Nir\\TriC\\HM_models_of_units\\1A6D\\Atemp_HM_"+units.charAt(d)+".pdb");
//				PutUnitInAnyTopPosition.alignByEquatorialDomains(pos1, 'I', alLeft, units.charAt(c));
//				PutUnitInAnyTopPosition.alignByEquatorialDomains(pos2, 'I', alRight, units.charAt(d));
//				alLeft.setChain(""+units.charAt(c));
//				alRight.setChain(""+units.charAt(d));
//				alLeft.add(alRight);
//				xlVec.checkConsistencyWithStructure(alLeft,""+units.charAt(c),""+units.charAt(d));
//				System.out.println();
//			}
//		}

		
//		// Analysis of SASA
//		String upperSeq = "AGZQHEBD";
//		String lowerSeq = "HEBDAGZQ";
//		AtomList fullComplex = PutUnitInAnyTopPosition.buildFullComplex(upperSeq, lowerSeq);
//		DSSP dssp = new DSSP("fullComplex.dssp");
//		xlVec.printSASAstats(dssp, fullComplex);
//		System.exit(0);
		
		
//		// Analysis on a single arrangement.
////		CrosslinkVectorTRiC xlVec = new CrosslinkVectorTRiC("TRIC_ATP_xlinks_result_lowerconfidence_sorted_NK.txt",0);
////		CrosslinkVectorTRiC xlVec = new CrosslinkVectorTRiC("Abersold_Yeast_intra_inter.txt",0);
//		CrosslinkVectorTRiC xlVec = new CrosslinkVectorTRiC("../Bovine_Aebersold/Abersold_Bovine_intra_inter.txt",0);
//		CrosslinkVectorTRiC xlVec1 = new CrosslinkVectorTRiC("../Bovine_Aebersold/NK_combined_moderate.txt",1);
//		System.out.println(xlVec.compareToOtherXLvec(xlVec1));
//		System.out.println(xlVec.filterAboveScore(25.0).size() + " " + xlVec.filterAboveScore(25.0).filterOutInterUnit().size() + " " +
//				xlVec.filterAboveScore(25.0).filterOutIntraUnit().size());
//		System.exit(0);
////		CrosslinkVectorTRiC xlVecAdd = new CrosslinkVectorTRiC("NK_combined_moderate.txt",1);
////		xlVec.addVecNonRedundant(xlVecAdd);
//		CrosslinkVectorTRiC xlVecInter = xlVec;//(CrosslinkVectorTRiC) xlVec.filterOutIntraUnit();
//		String upperSeq = "AGZQHEBD";
//		String lowerSeq = "HEBDAGZQ";
//		xlVecInter.setViolationCutoff(28.0);
//		AtomList fullComplex = PutUnitInAnyTopPosition.buildFullComplex(upperSeq, lowerSeq);
//		int[] results = xlVecInter.checkConsistencyWithStructure(fullComplex, "X", "X", true);
//		try {
//			fullComplex.filter(new AtomList.NonHydrogen()).print(new MeshiWriter("fullComplex.pdb"));
//		} catch (IOException e) {
//			e.printStackTrace();
//		}

		
//		// Analysis on the XL scoring according to Aebersold. This part was meant to see at what 
//		// rank (of their confidence) serious vioations start to occur. In order for it to work you need
//		// to establish a 'toTake' field and then to run the loop on crosslinks in: 'checkConsistencyWithStructure'.
//		// from 0 to 'toTake'.
//		CrosslinkVector xlVecInter = new CrosslinkVector("TRIC_xlinks_result_lowerconfidence.txt");
//		String upperSeq = "AGZQHEBD";
//		String lowerSeq = "HEBDAGZQ";
//		xlVecInter.setViolationCutoff(28.0);
//		for (int toTakeC=1 ; toTakeC<150 ; toTakeC++) {
//		xlVecInter.takeTo = toTakeC;
//		AtomList fullComplex = PutUnitInAnyTopPosition.buildFullComplex(upperSeq, lowerSeq);
//		int[] results = xlVecInter.checkConsistencyWithStructure(fullComplex, "X", "X", false);
//		System.out.println(toTakeC + " " + results[0] + " " + results[1] + " " + results[2]);
//		}

		
//		// This part calculates the violations/consistencies of all the 8! arrangements 
//		// against one Crosslink Vector.
//		// ----------------------------------------------------------------------------
//		CrosslinkVectorTRiC xlVec = (CrosslinkVectorTRiC) (new CrosslinkVectorTRiC("TRIC_ATP_xlinks_result_sorted_NK.txt",0)).filterOutRedundancies();
//		int firstPermute = Integer.parseInt(args[0]);
//		int lastPermute = Integer.parseInt(args[1]);
//		double violationCutoff = 27.55;
//		if (args.length>2) {
//			violationCutoff = Double.parseDouble(args[2]);
//		}
//		String[] outUpper = new String[8*(lastPermute-firstPermute+1)];
//		String[] outLower = new String[8*(lastPermute-firstPermute+1)];
//		String[] outScore = new String[8*(lastPermute-firstPermute+1)];
//		String units = "ABGDEHQZ";
//		CrosslinkVectorTRiC xlVecInter = (CrosslinkVectorTRiC) xlVec.filterOutIntraUnit();
//		xlVecInter.setViolationCutoff(violationCutoff);
//		String[] ringPerms = File2StringArray.f2a("ringPermute.txt");
//		int counter = 0;
//		for (int perm = firstPermute ; perm<=lastPermute ; perm++) {
//			for (int reg = 0 ; reg<8 ; reg++) {
//				String upperSeq = "";
//				String lowerSeq = "";
//				for (int c=0 ; c<8 ; c++) {
//					int permAtCup = Integer.parseInt(""+ringPerms[perm].charAt(c))-1;
//					int permAtCdown = Integer.parseInt(""+ringPerms[perm].charAt((reg + c) % 8))-1;
//					upperSeq += units.charAt(permAtCup);
//					lowerSeq += units.charAt(permAtCdown);
//				}
//				System.out.println("\n\n" + upperSeq +"\n" + lowerSeq);
//				AtomList fullComplex = PutUnitInAnyTopPosition.buildFullComplex(upperSeq, lowerSeq);
//				System.out.print("999999 ");
//				int[] result = xlVecInter.checkConsistencyWithStructure(fullComplex, "X", "X", false);
//				System.out.println();
//				System.out.println("Did perm:" + perm + "    Reg: " + reg);
//				outUpper[counter] = upperSeq;
//				outLower[counter] = lowerSeq;
//				outScore[counter] = result[0] + " " + result[1] + " " + result[2];
//				counter++;
//			}
//		}
//		// Writing the results to disk.
//		try{
//			BufferedWriter bwUp = new BufferedWriter(new FileWriter("xl_"+firstPermute+"_"+lastPermute+"_"+violationCutoff+"_upper.txt"));
//			BufferedWriter bwLow = new BufferedWriter(new FileWriter("xl_"+firstPermute+"_"+lastPermute+"_"+violationCutoff+"_lower.txt"));
//			BufferedWriter bwScore = new BufferedWriter(new FileWriter("xl_"+firstPermute+"_"+lastPermute+"_"+violationCutoff+"_score.txt"));
//			for (int c=0 ; c<outUpper.length ; c++) {
//				bwUp.write(outUpper[c] + "\n");
//				bwLow.write(outLower[c] + "\n");
//				bwScore.write(outScore[c] + "\n");
//			}
//			bwUp.close();
//			bwLow.close();
//			bwScore.close();
//		}
//		catch(Exception e) {
//		    throw new RuntimeException(e.getMessage());
//		}  		


//		// This part calculates the sum-violation score of all the 8! arrangements 
//		// against one Crosslink Vector.
//		// ----------------------------------------------------------------------------
//		CrosslinkVectorTRiC xlVec = (CrosslinkVectorTRiC) 
//				new CrosslinkVectorTRiC("Old_2010_set_Aebersold_Bovine.txt",0).filterOutRedundancies().filterOutIntraUnit();
//		System.out.println("Total number of XLs:" + xlVec.size());
//		int firstPermute = Integer.parseInt(args[0]);
//		int lastPermute = Integer.parseInt(args[1]);
//		double violationCutoff = 24.0;
//		if (args.length>2) {
//			violationCutoff = Double.parseDouble(args[2]);
//		}
//		String[] outUpper = new String[8*(lastPermute-firstPermute+1)];
//		String[] outLower = new String[8*(lastPermute-firstPermute+1)];
//		String[] outScore = new String[8*(lastPermute-firstPermute+1)];
//		String units = "ABGDEHQZ";
//		CrosslinkVectorTRiC xlVecInter = (CrosslinkVectorTRiC) xlVec.filterOutIntraUnit();
//		xlVecInter.setViolationCutoff(violationCutoff);
//		String[] ringPerms = File2StringArray.f2a("ringPermute.txt");
//		int counter = 0;
//		for (int perm = firstPermute ; perm<=lastPermute ; perm++) {
//			for (int reg = 0 ; reg<8 ; reg++) {
//				String upperSeq = "";
//				String lowerSeq = "";
//				for (int c=0 ; c<8 ; c++) {
//					int permAtCup = Integer.parseInt(""+ringPerms[perm].charAt(c))-1;
//					int permAtCdown = Integer.parseInt(""+ringPerms[perm].charAt((reg + c) % 8))-1;
//					upperSeq += units.charAt(permAtCup);
//					lowerSeq += units.charAt(permAtCdown);
//				}
//				System.out.println("\n\n" + upperSeq +"\n" + lowerSeq);
//				AtomList fullComplex = PutUnitInAnyTopPosition.buildFullComplex(upperSeq, lowerSeq);
//				System.out.print("999999 ");
//				double[] result = xlVecInter.checkConsistencyWithStructureSumOfViolations(fullComplex, "X", "X", false);
//				System.out.println();
//				System.out.println(result[0] + " " + result[1] + " " + result[2] + " " + result[3]);
//				System.out.println("Did perm:" + perm + "    Reg: " + reg);
//				outUpper[counter] = upperSeq;
//				outLower[counter] = lowerSeq;
//				outScore[counter] = result[0] + " " + result[1] + " " + result[2] + " " + result[3];
//				counter++;
//			}
//		}
//		// Writing the results to disk.
//		try{
//			BufferedWriter bwUp = new BufferedWriter(new FileWriter("xl_"+firstPermute+"_"+lastPermute+"_"+violationCutoff+"_upper.txt"));
//			BufferedWriter bwLow = new BufferedWriter(new FileWriter("xl_"+firstPermute+"_"+lastPermute+"_"+violationCutoff+"_lower.txt"));
//			BufferedWriter bwScore = new BufferedWriter(new FileWriter("xl_"+firstPermute+"_"+lastPermute+"_"+violationCutoff+"_score.txt"));
//			for (int c=0 ; c<outUpper.length ; c++) {
//				bwUp.write(outUpper[c] + "\n");
//				bwLow.write(outLower[c] + "\n");
//				bwScore.write(outScore[c] + "\n");
//			}
//			bwUp.close();
//			bwLow.close();
//			bwScore.close();
//		}
//		catch(Exception e) {
//		    throw new RuntimeException(e.getMessage());
//		}  		

		
		
//		// The Shuffling experiment 
//		CrosslinkVector xlVecTemp = new CrosslinkVector("TRIC_ATP_xlinks_result.txt");		
//		CrosslinkVector xlVec = new CrosslinkVector("TRIC_ATP_xlinks_result.txt");		
//		int firstPermute = Integer.parseInt(args[0]);
//		int lastPermute = Integer.parseInt(args[1]);
//		int permNum = Integer.parseInt(args[2]);
//
//		String[] outScore = new String[8*(lastPermute-firstPermute+1)];
//		String units = "ABGDEHQZ";
//		CrosslinkVector xlVecInterTemp = xlVecTemp.filterOutIntraUnit();
//		CrosslinkVector xlVecInter = xlVec.filterOutIntraUnit();
//		String[] ringPerms = File2StringArray.f2a("ringPermute.txt");
//
//		// The shuffling bit - starts
//		String[] randomPerm = File2StringArray.f2a("XL_Permutations.txt");
//		String thePerm = randomPerm[permNum];
//		StringTokenizer stPerm = new StringTokenizer(thePerm);
//		for (int c = 0 ; c<xlVecInter.size() ; c++) {
//			int take2from = Integer.parseInt(stPerm.nextToken())-1;
//			System.out.print(take2from + " ");
//			xlVecInter.get(c).setAbsPos2(xlVecInterTemp.get(take2from).absPos2());
//			xlVecInter.get(c).setProtName2(xlVecInterTemp.get(take2from).protName2());
//		}
//		System.out.println();
//		// The shuffling bit - ends
//		
//		int counter = 0;
//		for (int perm = firstPermute ; perm<=lastPermute ; perm++) {
//			for (int reg = 0 ; reg<8 ; reg++) {
//				String upperSeq = "";
//				String lowerSeq = "";
//				for (int c=0 ; c<8 ; c++) {
//					int permAtCup = Integer.parseInt(""+ringPerms[perm].charAt(c))-1;
//					int permAtCdown = Integer.parseInt(""+ringPerms[perm].charAt((reg + c) % 8))-1;
//					upperSeq += units.charAt(permAtCup);
//					lowerSeq += units.charAt(permAtCdown);
//				}
//				System.out.println("\n\n" + upperSeq +"\n" + lowerSeq);
//				AtomList fullComplex = PutUnitInAnyTopPosition.buildFullComplex(upperSeq, lowerSeq);
//				int[] result = xlVecInter.checkConsistencyWithStructure(fullComplex, "X", "X", false);
//				System.out.println("Did perm:" + perm + "    Reg: " + reg);
//				outScore[counter] = result[0] + " " + result[1] + " " + result[2];
//				counter++;
//			}
//		}
//		// Writing the results to disk.
//		try{
//			BufferedWriter bwScore = new BufferedWriter(new FileWriter("xl_"+firstPermute+"_"+lastPermute+"_"+permNum+"_score.txt"));
//			for (int c=0 ; c<outScore.length ; c++) {
//				bwScore.write(outScore[c] + "\n");
//			}
//			bwScore.close();
//		}
//		catch(Exception e) {
//		    throw new RuntimeException(e.getMessage());
//		}  				
		
		

//		// This part was meant to print a list of all close lysines per arrangement. 
//		// I did it at Michael's request. 
//		// This part does not really use any cross linking data.
//		int firstPermute = Integer.parseInt(args[0]);
//		int lastPermute = Integer.parseInt(args[1]);
//		double lysineCutoff = 30.0;
//		int baseLine = (firstPermute*8);
//		String upperSeqFilename = "upperSeq.txt";
//		String lowerSeqFilename = "lowerSeq.txt";
//		String closeLysFilename = "closeLys.txt";
//		try{
//			BufferedWriter bwCluseLys = new BufferedWriter(new FileWriter("xl_"+firstPermute+"_"+lastPermute+"_"+closeLysFilename));
//			String[] outUpper = new String[8*(lastPermute-firstPermute+1)];
//			String[] outLower = new String[8*(lastPermute-firstPermute+1)];
//			String units = "ABGDEHQZ";
//			String[] ringPerms = File2StringArray.f2a("ringPermute.txt");
//			int counter = 0;
//			for (int perm = firstPermute ; perm<=lastPermute ; perm++) {
//				for (int reg = 0 ; reg<8 ; reg++) {
//					String upperSeq = "";
//					String lowerSeq = "";
//					String lysDis = "";
//
//					for (int c=0 ; c<8 ; c++) {
//						int permAtCup = Integer.parseInt(""+ringPerms[perm].charAt(c))-1;
//						int permAtCdown = Integer.parseInt(""+ringPerms[perm].charAt((reg + c) % 8))-1;
//						upperSeq += units.charAt(permAtCup);
//						lowerSeq += units.charAt(permAtCdown);
//					}
//					// making the 'lowerSeq' readable in our common format
//					String semiInverseLowerSeq = "" + lowerSeq.charAt(0);
//					for (int tmpc=7 ; tmpc>0 ; tmpc--) {
//						semiInverseLowerSeq += lowerSeq.charAt(tmpc);
//					}
//					System.out.println("\n\n" + upperSeq +"\n" + semiInverseLowerSeq);
//					AtomList fullComplex = PutUnitInAnyTopPosition.buildFullComplex(upperSeq, lowerSeq);
//					AtomList justKCB = fullComplex.filter(new AtomList.KCB_Filter());
//					for (int c1=0 ; c1<justKCB.size() ; c1++) {
//						for (int c2=c1+1 ; c2<justKCB.size() ; c2++) {
//							if (justKCB.atomAt(c1).distanceFrom(justKCB.atomAt(c2)) < lysineCutoff) {
//								lysDis += ((baseLine+counter) + " " + justKCB.atomAt(c1).chain() + " " + justKCB.atomAt(c1).residueNumber() + " " +  
//										justKCB.atomAt(c2).chain() + " " + justKCB.atomAt(c2).residueNumber() + " " + ((int) justKCB.atomAt(c1).distanceFrom(justKCB.atomAt(c2))) + "\n");
//							}
//						}					
//					}
//					bwCluseLys.write(lysDis);				
//					System.out.println("Did perm:" + perm + "    Reg: " + reg);
//					outUpper[counter] = upperSeq;
//					outLower[counter] = semiInverseLowerSeq;
//					counter++;
//				}
//			}
//			// Writing the arrangement to disk.
//			BufferedWriter bwUp = new BufferedWriter(new FileWriter("xl_"+firstPermute+"_"+lastPermute+"_"+upperSeqFilename));
//			BufferedWriter bwLow = new BufferedWriter(new FileWriter("xl_"+firstPermute+"_"+lastPermute+"_"+lowerSeqFilename));
//			for (int c=0 ; c<outUpper.length ; c++) {
//				bwUp.write(outUpper[c] + "\n");
//				bwLow.write(outLower[c] + "\n");
//			}
//			bwUp.close();
//			bwLow.close();
//			bwCluseLys.close();
//		}
//		catch(Exception e) {
//			throw new RuntimeException(e.getMessage());
//		}  		

		
//		// Printing a XL list with stuctured/unstructured information.
//		CrosslinkVectorTRiC xlVec = (CrosslinkVectorTRiC) (new CrosslinkVectorTRiC("TRIC_ATP_xlinks_result_sorted_NK.txt",0));
////		CrosslinkVectorTRiC xlVecNir = (CrosslinkVectorTRiC) (new CrosslinkVectorTRiC("NK_combined_moderate.txt",1));
////		xlVec.addVecNonRedundant(xlVecNir);
//		String upperSeq = "AGZQHEBD";
//		String lowerSeq = "HEBDAGZQ";
//		AtomList fullComplex = PutUnitInAnyTopPosition.buildFullComplex(upperSeq, lowerSeq);
//		System.out.print(xlVec.printListWithStructuralInfo(fullComplex));


//		// This part calculates the violations/consistencies of all the 7! ring arrangements 
//		// against one Crosslink Vector. Don't forget to disable the bottom ring in 'identicalChains'.  
//		// -------------------------------------------------------------------------------------------
//		CrosslinkVectorTRiC xlVec = new CrosslinkVectorTRiC("TRIC_ATP_xlinks_result_sorted_NK.txt",0);
//		CrosslinkVectorTRiC xlVecAdd = new CrosslinkVectorTRiC("NK_combined_moderate.txt",1);
//		xlVec.addVecNonRedundant(xlVecAdd);		
//		int firstPermute = Integer.parseInt(args[0]);
//		int lastPermute = Integer.parseInt(args[1]);
//		double violationCutoff = 27.55;
//		if (args.length>2) {
//			violationCutoff = Double.parseDouble(args[2]);
//		}
//		String[] outUpper = new String[(lastPermute-firstPermute+1)];
//		String[] outScore = new String[(lastPermute-firstPermute+1)];
//		String units = "ABGDEHQZ";
//		CrosslinkVectorTRiC xlVecInter = (CrosslinkVectorTRiC) xlVec.filterOutIntraUnit();
//		xlVecInter.setViolationCutoff(violationCutoff);
//		String[] ringPerms = File2StringArray.f2a("ringPermute.txt");
//		int counter = 0;
//		for (int perm = firstPermute ; perm<=lastPermute ; perm++) {
//			String upperSeq = "";
//			for (int c=0 ; c<8 ; c++) {
//				int permAtCup = Integer.parseInt(""+ringPerms[perm].charAt(c))-1;
//				upperSeq += units.charAt(permAtCup);
//			}
//			System.out.println("\n\n" + upperSeq);
//			AtomList fullComplex = PutUnitInAnyTopPosition.buildFullComplex(upperSeq, upperSeq);
//			int[] result = xlVecInter.checkConsistencyWithStructure(fullComplex, "X", "X", false);
//			System.out.println("Did perm:" + perm);
//			outUpper[counter] = upperSeq;
//			outScore[counter] = result[0] + " " + result[1] + " " + result[2];
//			counter++;
//		}
//		// Writing the results to disk.
//		try{
//			BufferedWriter bwUp = new BufferedWriter(new FileWriter("xl_"+firstPermute+"_"+lastPermute+"_"+violationCutoff+"_upper.txt"));
//			BufferedWriter bwLow = new BufferedWriter(new FileWriter("xl_"+firstPermute+"_"+lastPermute+"_"+violationCutoff+"_lower.txt"));
//			BufferedWriter bwScore = new BufferedWriter(new FileWriter("xl_"+firstPermute+"_"+lastPermute+"_"+violationCutoff+"_score.txt"));
//			for (int c=0 ; c<outUpper.length ; c++) {
//				bwUp.write(outUpper[c] + "\n");
//				bwScore.write(outScore[c] + "\n");
//			}
//			bwUp.close();
//			bwLow.close();
//			bwScore.close();
//		}
//		catch(Exception e) {
//		    throw new RuntimeException(e.getMessage());
//		}  		
		
//		// This part will print the distance of each cross-link on the model. 
////		CrosslinkVectorTRiC xlVec = new CrosslinkVectorTRiC("NK_combined_moderate.txt",1);
//		CrosslinkVectorTRiC xlVec = new CrosslinkVectorTRiC("Abersold_Yeast_intra_inter.txt",0);
//		String upperSeq = "AGZQHEBD";
//		String lowerSeq = "HEBDAGZQ";
//		AtomList fullComplex = PutUnitInAnyTopPositionYeast.buildFullComplex(upperSeq, lowerSeq);
//		for (Crosslink xl : xlVec) {
//			String str = "N/A";
//			if ((fullComplex.findAtomInListReturningAtom("CA", xl.protName1(), xl.absPos1()) != null) &
//				(fullComplex.findAtomInListReturningAtom("CA", ""+xl.protName2().charAt(0), xl.absPos2()) != null)) {
//				double sameRingDis = fullComplex.findAtomInListReturningAtom("CA", xl.protName1(), xl.absPos1()).distanceFrom(
//						fullComplex.findAtomInListReturningAtom("CA", ""+xl.protName2().charAt(0), xl.absPos2()) );
//				double acrossRingDis = fullComplex.findAtomInListReturningAtom("CA", xl.protName1(), xl.absPos1()).distanceFrom(
//						fullComplex.findAtomInListReturningAtom("CA", ""+PutUnitInAnyTopPositionYeast.matchChainLetterOnBottomRing(xl.protName2().charAt(0)), xl.absPos2()) );
//				if (Math.min(sameRingDis, acrossRingDis)<0.01) {
//					str = ""+ Math.min(sameRingDis, acrossRingDis) + " " + Math.max(sameRingDis, acrossRingDis) ;
//				}
//				else {
//					str = ""+Math.min(sameRingDis, acrossRingDis);
//				}
//				System.out.println(str);
//			}
//		}
		

//		// This part compare new results from Feb 2015 with Aebersold and Leitner paper. 
//		CrosslinkVector xlVec = new CrosslinkVectorTRiC("C:\\Users\\Nir\\CCT_HUJI\\2015_02_18_analysis\\Report_2015_02_18_01_CCT-XL_new_above055.txt",1);
//		CrosslinkVector other = new CrosslinkVectorTRiC("C:\\Users\\Nir\\CCT_HUJI\\2015_02_18_analysis\\Yeast_APO_cross-links_Leitner.txt",0);
//		CrosslinkVector xlVec = new CrosslinkVectorTRiC("C:\\Users\\Nir\\CCT_HUJI\\2015_02_18_analysis\\Report_2015_02_18_04_CCT-XL_old_above055.txt",1);
//		String upperSeq = "AGZQHEBD";
//		String lowerSeq = "HEBDAGZQ";
//		AtomList fullComplex = PutUnitInAnyTopPositionYeast.buildFullComplex(upperSeq, lowerSeq);
//		System.out.print(xlVec.printListWithStructuralInfo(fullComplex));
//
//		CrosslinkVector xlVec1 = new CrosslinkVectorTRiC("C:\\Users\\Nir\\CCT_HUJI\\2015_02_18_analysis\\Report_2015_02_18_01_CCT-XL_new_above055.txt",1);
//		CrosslinkVector xlVec2 = new CrosslinkVectorTRiC("C:\\Users\\Nir\\CCT_HUJI\\2015_02_18_analysis\\Report_2015_02_22_15_CCT-XL_new_above055.txt",1);
//		CrosslinkVector xlVec3 = new CrosslinkVectorTRiC("C:\\Users\\Nir\\CCT_HUJI\\2015_02_18_analysis\\Report_2015_02_22_17_CCT-XL_old_above055.txt",1);
//		xlVec.addVecNonRedundant(xlVec1);
//		xlVec.addVecNonRedundant(xlVec2);
//		xlVec.addVecNonRedundant(xlVec3);
//		CrosslinkVector xlVec = new CrosslinkVectorTRiC("C:\\Users\\Nir\\CCT_HUJI\\2015_02_18_analysis\\Report_2015_02_18_04_CCT-XL_old_decoy_1_above055.txt",1);
//		CrosslinkVector Rawother = new CrosslinkVectorTRiC("C:\\Users\\Nir\\CCT_HUJI\\2015_02_18_analysis\\Yeast_APO_cross-links_Leitner.txt",0);
//		CrosslinkVector other = Rawother.filterOutRedundancies().filterOutIntraUnit();
//		System.out.println("Size us:" + xlVec.size() + " Inter:" + xlVec.filterOutIntraUnit().size() + "\nSize Aeber:" + other.size() + " Inter:" + other.filterOutIntraUnit().size() );
//		System.out.println(xlVec.compareToOtherXLvec(other));
		
		
		
		
//		// This part will print the histogram of all possible Lysine-Lysine distances. 
//		int maxDis = 300;
//		double[] occurred = new double[maxDis];
//		String upperSeq = "AGZQHEBD";
//		String lowerSeq = "HEBDAGZQ";
//		AtomList fullComplex = PutUnitInAnyTopPosition.buildFullComplex(upperSeq, lowerSeq);
//		for (int c1=0 ; c1<fullComplex.size()/2 ; c1++) {
//			if (c1/1000 == c1/1000.0){
//				System.out.println(c1);
//			}
//			if (fullComplex.atomAt(c1).name().equals("CA") && fullComplex.atomAt(c1).residueName().equals("LYS")) {
//				for (int c2=c1+1 ; c2<fullComplex.size()/2 ; c2++) {
//					if (fullComplex.atomAt(c2).name().equals("CA") && fullComplex.atomAt(c2).residueName().equals("LYS")) {
//						double dis1 = fullComplex.atomAt(c2).distanceFrom(fullComplex.atomAt(c1));
//						Atom atom2 = fullComplex.findAtomInListReturningAtom("CA", 
//								PutUnitInAnyTopPosition.matchChainLetterOnBottomRing(fullComplex.atomAt(c2).chain().charAt(0))+"" ,
//								fullComplex.atomAt(c2).residueNumber());
//						double dis2 = atom2.distanceFrom(fullComplex.atomAt(c1));
//						double dis = Math.min(dis1, dis2);
//						occurred[(int) (dis/5.0)]++;
//					}
//				}
//			}
//		}
//		System.out.println("****************************************");
//		for (int c1=0 ; c1<maxDis ; c1+=5) {
//			System.out.println(c1+" " + occurred[c1/5]);
//		}

		
//		// This part will put crosslinker atoms in place. 
//		String upperSeq = "AGZQHEBD";
//		String lowerSeq = "HEBDAGZQ";
//		AtomList fullComplex = PutUnitInAnyTopPosition.buildFullComplex(upperSeq, lowerSeq);
//		PuttingXLinPDB xlInPDB = new PuttingXLinPDB(fullComplex);
//		CrosslinkVectorTRiC xlVec = (CrosslinkVectorTRiC) (new CrosslinkVectorTRiC("NK_combined_moderate_for_Fig4.txt",1)).filterOutIntraUnit();
////		xlInPDB.insideOutsideAllvec(xlVec);
//		for (int c=0 ; c<xlVec.size() ; c+=2) {
//			System.out.println(xlVec.elementAt(c));
//			xlInPDB.refineXLatoms(xlVec.elementAt(c) , c);
//			switch (xlVec.elementAt(c).protName1().charAt(0)) {
//			   case 'A':  xlVec.elementAt(c).setProtName1("I"); break;
//			   case 'B':  xlVec.elementAt(c).setProtName1("J"); break;
//			   case 'G':  xlVec.elementAt(c).setProtName1("K"); break;
//			   case 'D':  xlVec.elementAt(c).setProtName1("L"); break;
//			   case 'E':  xlVec.elementAt(c).setProtName1("M"); break;
//			   case 'H':  xlVec.elementAt(c).setProtName1("N"); break;
//			   case 'Q':  xlVec.elementAt(c).setProtName1("O"); break;
//			   case 'Z':  xlVec.elementAt(c).setProtName1("P"); break;
//			   case 'J':  xlVec.elementAt(c).setProtName2("B"); break;
//			   case 'L':  xlVec.elementAt(c).setProtName2("D"); break;
//			   case 'M':  xlVec.elementAt(c).setProtName1("E"); break;
//			   case 'P':  xlVec.elementAt(c).setProtName1("Z"); break;
//			}
//			switch (xlVec.elementAt(c).protName2().charAt(0)) {
//			   case 'A':  xlVec.elementAt(c).setProtName2("I"); break;
//			   case 'B':  xlVec.elementAt(c).setProtName2("J"); break;
//			   case 'G':  xlVec.elementAt(c).setProtName2("K"); break;
//			   case 'D':  xlVec.elementAt(c).setProtName2("L"); break;
//			   case 'E':  xlVec.elementAt(c).setProtName2("M"); break;
//			   case 'H':  xlVec.elementAt(c).setProtName2("N"); break;
//			   case 'Q':  xlVec.elementAt(c).setProtName2("O"); break;
//			   case 'Z':  xlVec.elementAt(c).setProtName2("P"); break;
//			   case 'J':  xlVec.elementAt(c).setProtName2("B"); break;
//			   case 'L':  xlVec.elementAt(c).setProtName2("D"); break;
//			   case 'M':  xlVec.elementAt(c).setProtName1("E"); break;
//			   case 'P':  xlVec.elementAt(c).setProtName1("Z"); break;
//			}
//			System.out.println(xlVec.elementAt(c));
//			xlInPDB.refineXLatoms(xlVec.elementAt(c) , c+1);
//		}		
	} // Of main()
	
	
}
