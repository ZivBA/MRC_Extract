 package meshi.applications.TriC.yeast;

import java.io.IOException;

import meshi.applications.TriC.TricYeastAlignment;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.file.MeshiWriter;
import meshi.util.overlap.Overlap;

public class OverlapSubunits {

	// This "main" will overlap the right unit in each of the eight pairs on the equatorial domain of 1Q3R	
	public static void main(String[] args) {
		// This code will overlap the ring pairs AG, GZ, etc... with the 1Q3R chain A.
		int[][] alignBy = {{79,90},
				{123,140},
				{501,514}}; // In 1Q3R frame
		String refFile = "C:/Users/Nir/TRiC/Crystallography/Intra_ring_interface/1Q3R.pdb";
//		int[][] alignBy = {{101,112},
//				{145,161},
//				{534,547}}; // In 1Q3R frame
//		String refFile = "C:/Users/Nir/TRiC/Crystallography/Intra_ring_interface/Manual_optimization_on_LJcap_15/refine_16-half.ENCAD.refined_interface.pdb";
		String ringFile = "C:/Users/Nir/TRiC/Crystallography/Intra_ring_interface/Manual_optimization_on_LJcap_15/refine_16-half.ENCAD.refined_interface.MichaelChain.pdb";
		String firstLetterChain = "AGZQHEBD";
		String secondLetterChain = "GZQHEBDA";
		String firstLetterGene = "AGZQHEBD";
		String secondLetterGene = "GZQHEBDA";
		AtomList ref = new AtomList(refFile);
		ref.chainFilter("A");
//		ref.chainFilter("E");
		int commonAtomNum = 0;
		for (int rangeC=0 ; rangeC<alignBy.length ; rangeC++) {
			for (int resC=alignBy[rangeC][0] ; resC<=alignBy[rangeC][1] ; resC++) {
				commonAtomNum++;
			}			
		}
		int[] refResNums = new int[commonAtomNum];
		int counter = 0;
		for (int rangeC=0 ; rangeC<alignBy.length ; rangeC++) {
			for (int resC=alignBy[rangeC][0] ; resC<=alignBy[rangeC][1] ; resC++) {
				refResNums[counter] = resC;
				counter++;
			}			
		}
		TricYeastAlignment alignment = new TricYeastAlignment();
		for (int c=0 ; c<8 ; c++) {
			Atom.resetNumberOfAtoms();
			AtomList ring = new AtomList(ringFile);
			AtomList prot1 = ring.chainFilter(firstLetterChain.charAt(c)+"");
			AtomList prot2 = ring.chainFilter(secondLetterChain.charAt(c)+"");
			prot1.setChain(""+firstLetterGene.charAt(c));
			prot2.setChain(""+secondLetterGene.charAt(c));
			AtomList pair = new AtomList();
			pair.add(filterDomainAccordingTo1Q3R(prot2, secondLetterGene.charAt(c), 'E'));
			pair.add(filterDomainAccordingTo1Q3R(prot1, firstLetterGene.charAt(c), 'E'));
			int[] moveResNums = new int[commonAtomNum];
			counter = 0;
			for (int rangeC=0 ; rangeC<alignBy.length ; rangeC++) {
				for (int resC=alignBy[rangeC][0] ; resC<=alignBy[rangeC][1] ; resC++) {
					moveResNums[counter] = alignment.getNewResNum('K', resC, secondLetterGene.charAt(c));
//					moveResNums[counter] = alignment.getNewResNum('E', resC, secondLetterGene.charAt(c));
//					System.out.print(alignment.getNewResNum('E', resC, secondLetterGene.charAt(c))+",");
					counter++;
				}
//				System.out.println();
			}
			double rms = Overlap.fitProteins(ref, refResNums, pair, moveResNums);
			System.out.println("The subunit " + secondLetterGene.charAt(c) + " overlap is: " + ((int) (rms*100))/100.0);
			try {
				pair.print(new MeshiWriter("C:/Users/Nir/TRiC/Crystallography/Figures/Middle_domain_position/PAIR_" + firstLetterGene.charAt(c) + 
						secondLetterGene.charAt(c)+ ".pdb"));
			} catch (IOException e) {
				throw new RuntimeException("Failed to write chains");
			}
		}
	}
	

// This "main" will overlap the left unit in each of the eight pairs on the equatorial domain of 1Q3R	
//	public static void main(String[] args) {
//		// This code will overlap the ring pairs AG, GZ, etc... with the 1Q3R chain A.
//		int[][] alignBy = {{24,39},
//				{79,90},
//				{99,115},
//				{123,140},
//				{415,424},
//				{435,453},
//				{461,471},
//				{501,514}}; // In 1Q3R frame
//		String refFile = "1Q3R.pdb";
//		String ringFile = "refine_16.pdb";
//		String firstLetterChain = "AGDHCBFE";
//		String secondLetterChain = "GDHCBFEA";
//		String firstLetterGene = "AGZQHEBD";
//		String secondLetterGene = "GZQHEBDA";
//		AtomList ref = new AtomList(refFile);
//		ref.chainFilter("A");
//		int commonAtomNum = 0;
//		for (int rangeC=0 ; rangeC<alignBy.length ; rangeC++) {
//			for (int resC=alignBy[rangeC][0] ; resC<=alignBy[rangeC][1] ; resC++) {
//				commonAtomNum++;
//			}			
//		}
//		int[] refResNums = new int[commonAtomNum];
//		int counter = 0;
//		for (int rangeC=0 ; rangeC<alignBy.length ; rangeC++) {
//			for (int resC=alignBy[rangeC][0] ; resC<=alignBy[rangeC][1] ; resC++) {
//				refResNums[counter] = resC;
//				counter++;
//			}			
//		}
//		TricYeastAlignment alignment = new TricYeastAlignment();
//		for (int c=0 ; c<8 ; c++) {
//			Atom.resetNumberOfAtoms();
//			AtomList ring = new AtomList(ringFile);
//			AtomList prot1 = ring.chainFilter(firstLetterChain.charAt(c)+"");
//			AtomList prot2 = ring.chainFilter(secondLetterChain.charAt(c)+"");
//			prot1.setChain(""+firstLetterGene.charAt(c));
//			prot2.setChain(""+secondLetterGene.charAt(c));
//			AtomList pair = new AtomList();
//			pair.add(filterDomainAccordingTo1Q3R(prot1, firstLetterGene.charAt(c), 'E'));
//			pair.add(filterDomainAccordingTo1Q3R(prot2, secondLetterGene.charAt(c), 'E'));
//			int[] moveResNums = new int[commonAtomNum];
//			counter = 0;
//			for (int rangeC=0 ; rangeC<alignBy.length ; rangeC++) {
//				for (int resC=alignBy[rangeC][0] ; resC<=alignBy[rangeC][1] ; resC++) {
//					moveResNums[counter] = alignment.getNewResNum('K', resC, firstLetterGene.charAt(c));
//					counter++;
//				}			
//			}
//			double rms = Overlap.fitProteins(ref, refResNums, pair, moveResNums);
//			System.out.println("The subunit " + firstLetterGene.charAt(c) + " overlap is: " + ((int) (rms*100))/100.0);
//			try {
//				pair.print(new MeshiWriter("PAIR_" + firstLetterGene.charAt(c) + 
//						secondLetterGene.charAt(c)+ ".pdb"));
//			} catch (IOException e) {
//				throw new RuntimeException("Failed to write chains");
//			}
//		}
//	}

//	protected static AtomList filterDomainAccordingTo1Q3R(AtomList fullList, char subunit , char domain) {
//		int[][] parsingQ3R; 
//		TricYeastAlignment alignments = new TricYeastAlignment();
//		switch (domain) {
//    	case 'E': // Equatorial 
//    		int[][] newParse = {{1,148} , {405,535}};
//    		parsingQ3R = newParse;
//    		break;
//    	case 'M': // Equatorial 
//    		int[][] newParse1 = {{149,217} , {369,404}};
//    		parsingQ3R = newParse1;
//    		break;
//    	case 'A': // Apical 
//    		int[][] newParse2 = {{218,368}};
//    		parsingQ3R = newParse2;
//    		break;
//    	default:
//    		throw new RuntimeException("Invalid domain letter {E,M,A}");
//		}
//		
//		int[][] parsing = new int[parsingQ3R.length][2];
//		for (int c=0 ; c<parsingQ3R.length ; c++) {
//			parsing[c][0] = alignments.getNewResNum('K', parsingQ3R[c][0], subunit);
//			System.out.print(parsing[c][0] + ",");
//			parsing[c][1] = alignments.getNewResNum('K', parsingQ3R[c][1], subunit);
//			System.out.print(parsing[c][1] + ",");
//		}	
//		System.out.println();
//		if (domain=='E') {
//			parsing[0][0] = 1; // To include N-terminal
//			//parsing[1][1] = 1000; // To include ATP 
//		}
//		return filterDomain(fullList, parsing);
//	}
	
//  This is for showing the equatorial and middle together
	protected static AtomList filterDomainAccordingTo1Q3R(AtomList fullList, char subunit , char domain) {
		int[][] parsingQ3R; 
		TricYeastAlignment alignments = new TricYeastAlignment();
		switch (domain) {
    	case 'E': // Equatorial 
    		int[][] newParse = {{1,217} , {369,535}};
    		parsingQ3R = newParse;
    		break;
    	case 'A': // Apical 
    		int[][] newParse2 = {{218,368}};
    		parsingQ3R = newParse2;
    		break;
    	default:
    		throw new RuntimeException("Invalid domain letter {E,A}");
		}
		
		int[][] parsing = new int[parsingQ3R.length][2];
		for (int c=0 ; c<parsingQ3R.length ; c++) {
			parsing[c][0] = alignments.getNewResNum('K', parsingQ3R[c][0], subunit);
			System.out.print(parsing[c][0] + ",");
			parsing[c][1] = alignments.getNewResNum('K', parsingQ3R[c][1], subunit);
			System.out.print(parsing[c][1] + ",");
		}	
		System.out.println();
		if (domain=='E') {
			parsing[0][0] = 1; // To include N-terminal
			//parsing[1][1] = 1000; // To include ATP 
		}
		return filterDomain(fullList, parsing);
	}
	
	protected static AtomList filterDomain(AtomList fullList, int[][] parsing) {
		AtomList newList = new AtomList();
		for (int c=0 ; c<fullList.size() ; c++) {
			int resNum = fullList.atomAt(c).residueNumber();
			for (int segID=0 ; segID<parsing.length ; segID++ ) {
				if ((resNum>=parsing[segID][0]) & (resNum<=parsing[segID][1])) {
					newList.add(fullList.atomAt(c));
				}
			}
		}
		return newList;
	}
	
	
}
