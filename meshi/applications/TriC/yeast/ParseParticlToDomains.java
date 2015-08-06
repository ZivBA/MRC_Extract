package meshi.applications.TriC.yeast;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import meshi.applications.TriC.TricYeastAlignment;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.file.File2StringArray;
import meshi.util.file.MeshiWriter;

public class ParseParticlToDomains {

	public ParseParticlToDomains() {}

	/**
	 * Parameter 1 - the file name
	 * Parameter 2 - Prefix of output
	 */
	public static void main(String[] args) {	
		String[] pairs = {"AG","GZ","ZQ","QH","HE","EB","BD","DA","AB"};
		String outputPrefix = "C:\\Users\\Nir\\TRiC\\Crystallography\\Intra_ring_interface_of_refine_19i\\Chimera_analysis\\";
		for (int unitC = 0; unitC<pairs.length ; unitC++) {
			AtomList list;
			Atom.resetNumberOfAtoms();
			char unitName1;
			char unitName2;
			if (unitC<8) {
				list = new AtomList("C:\\Users\\Nir\\TRiC\\Crystallography\\refine_19i-pdb-hkl\\refine_19i_chains_by_genes.pdb");
				unitName1 = pairs[unitC].charAt(0);
				unitName2 = pairs[unitC].charAt(1);
			}
			else {
				list = new AtomList("C:\\Users\\Nir\\TRiC\\Crystallography\\Intra_ring_interface\\1Q3R.pdb");
				unitName1 = 'K';
				unitName2 = 'K';				
			}
			// Doing single subunit 
			AtomList chainToWrite = filterDomainAccordingTo1Q3R(list.chainFilter(pairs[unitC].charAt(0)+""),unitName1,'E');
//			System.out.println(chainToWrite.CAFilter().size());
			try {
				chainToWrite.print(new MeshiWriter(outputPrefix+pairs[unitC]+"_LEFT_EQ.pdb"));
			} catch (IOException e) {
				throw new RuntimeException("Failed to write chain");
			}
			chainToWrite = filterDomainAccordingTo1Q3R(list.chainFilter(pairs[unitC].charAt(0)+""),unitName1,'M');
//			System.out.println(chainToWrite.CAFilter().size());
			try {
				chainToWrite.print(new MeshiWriter(outputPrefix+pairs[unitC]+"_LEFT_MI.pdb"));
			} catch (IOException e) {
				throw new RuntimeException("Failed to write chain");
			}
			chainToWrite = filterDomainAccordingTo1Q3R(list.chainFilter(pairs[unitC].charAt(0)+""),unitName1,'A');
//			System.out.println(chainToWrite.CAFilter().size());
			try {
				chainToWrite.print(new MeshiWriter(outputPrefix+pairs[unitC]+"_LEFT_AP.pdb"));
			} catch (IOException e) {
				throw new RuntimeException("Failed to write chain");
			}
			chainToWrite = filterDomainAccordingTo1Q3R(list.chainFilter(pairs[unitC].charAt(1)+""),unitName2,'E');
			System.out.println(chainToWrite.atomAt(0).residueNumber()-1);
//			System.out.println(chainToWrite.CAFilter().size());
			try {
				chainToWrite.print(new MeshiWriter(outputPrefix+pairs[unitC]+"_RIGHT_EQ.pdb"));
			} catch (IOException e) {
				throw new RuntimeException("Failed to write chain");
			}
			chainToWrite = filterDomainAccordingTo1Q3R(list.chainFilter(pairs[unitC].charAt(1)+""),unitName2,'M');
//			System.out.println(chainToWrite.CAFilter().size());
			try {
				chainToWrite.print(new MeshiWriter(outputPrefix+pairs[unitC]+"_RIGHT_MI.pdb"));
			} catch (IOException e) {
				throw new RuntimeException("Failed to write chain");
			}
			chainToWrite = filterDomainAccordingTo1Q3R(list.chainFilter(pairs[unitC].charAt(1)+""),unitName2,'A');
//			System.out.println(chainToWrite.CAFilter().size());
			try {
				chainToWrite.print(new MeshiWriter(outputPrefix+pairs[unitC]+"_RIGHT_AP.pdb"));
			} catch (IOException e) {
				throw new RuntimeException("Failed to write chain");
			}
			// Doing pairs
			chainToWrite = filterDomainAccordingTo1Q3R(list.chainFilter(pairs[unitC].charAt(0)+""),unitName1,'E');
			chainToWrite.add(filterDomainAccordingTo1Q3R(list.chainFilter(pairs[unitC].charAt(1)+""),unitName2,'E'));
//			System.out.println(chainToWrite.CAFilter().size());
			try {
				chainToWrite.print(new MeshiWriter(outputPrefix+pairs[unitC]+"_EQ_EQ.pdb"));
			} catch (IOException e) {
				throw new RuntimeException("Failed to write chain");
			}
			chainToWrite = filterDomainAccordingTo1Q3R(list.chainFilter(pairs[unitC].charAt(0)+""),unitName1,'M');
			chainToWrite.add(filterDomainAccordingTo1Q3R(list.chainFilter(pairs[unitC].charAt(1)+""),unitName2,'E'));
//			System.out.println(chainToWrite.CAFilter().size());
			try {
				chainToWrite.print(new MeshiWriter(outputPrefix+pairs[unitC]+"_MI_EQ.pdb"));
			} catch (IOException e) {
				throw new RuntimeException("Failed to write chain");
			}
			chainToWrite = filterDomainAccordingTo1Q3R(list.chainFilter(pairs[unitC].charAt(0)+""),unitName1,'A');
			chainToWrite.add(filterDomainAccordingTo1Q3R(list.chainFilter(pairs[unitC].charAt(1)+""),unitName2,'A'));
//			System.out.println(chainToWrite.CAFilter().size());
			try {
				chainToWrite.print(new MeshiWriter(outputPrefix+pairs[unitC]+"_AP_AP.pdb"));
			} catch (IOException e) {
				throw new RuntimeException("Failed to write chain");
			}
		}

		
//		String units = "AGZQHEBD";
//		String chains= "AGDHCBFE";
//		AtomList list = new AtomList(args[0]);
//		Atom.resetNumberOfAtoms();
//		for (int unitC = 0; unitC<units.length() ; unitC++) {
//			// Doing Equatorial
//			AtomList chainToWrite = filterDomainAccordingTo1Q3R(list.chainFilter(chains.charAt(unitC)+""),units.charAt(unitC),'E').duplicate();
//			System.out.println(chainToWrite.CAFilter().size());
//			chainToWrite.setChain(""+units.charAt(unitC));
//			try {
//				chainToWrite.print(new MeshiWriter(args[1]+"_"+units.charAt(unitC)+"_EQ.pdb"));
//			} catch (IOException e) {
//				throw new RuntimeException("Failed to write chain");
//			}
//			// Doing Apical
//			chainToWrite = filterDomainAccordingTo1Q3R(list.chainFilter(chains.charAt(unitC)+""),units.charAt(unitC),'A').duplicate();
//			System.out.println(chainToWrite.CAFilter().size());
//			chainToWrite.setChain(""+units.charAt(unitC));
//			try {
//				chainToWrite.print(new MeshiWriter(args[1]+"_"+units.charAt(unitC)+"_AP.pdb"));
//			} catch (IOException e) {
//				throw new RuntimeException("Failed to write chain");
//			}
//		}
		
//		String[] pairs = {"AG","GZ","ZQ","QH","HE","EB","BD","DA"};
//		for (int pc=0 ; pc<pairs.length ; pc++) {
//			char unit1 = pairs[pc].charAt(0);
//			char unit2 = pairs[pc].charAt(1);
//			String[] p1 = File2StringArray.f2a("domain_SCWRL_"+unit1+"_EQ.pdb");
//			String[] p2 = File2StringArray.f2a("domain_SCWRL_"+unit2+"_EQ.pdb");
//			try {
//				BufferedWriter bw = new BufferedWriter(new FileWriter("pair_SCWRL_"+pairs[pc]+"_EQ.pdb"));
//				for (int c=0; c<p1.length ; c++) {
//					bw.write(p1[c] + "\n");
//				}
//				bw.write("TER\n");
//				for (int c=0; c<p2.length ; c++) {
//					bw.write(p2[c] + "\n");
//				}
//				bw.write("TER\nEND\n");
//				bw.close();
//			} catch (IOException e) {
//				throw new RuntimeException("Failed to write chain");
//			}
//			p1 = File2StringArray.f2a("domain_SCWRL_"+unit1+"_AP.pdb");
//			p2 = File2StringArray.f2a("domain_SCWRL_"+unit2+"_AP.pdb");
//			try {
//				BufferedWriter bw = new BufferedWriter(new FileWriter("pair_SCWRL_"+pairs[pc]+"_AP.pdb"));
//				for (int c=0; c<p1.length ; c++) {
//					bw.write(p1[c] + "\n");
//				}
//				bw.write("TER\n");
//				for (int c=0; c<p2.length ; c++) {
//					bw.write(p2[c] + "\n");
//				}
//				bw.write("TER\nEND\n");
//				bw.close();
//			} catch (IOException e) {
//				throw new RuntimeException("Failed to write chain");
//			}
//		}

	}

	
	protected static AtomList filterDomainAccordingTo1Q3R(AtomList fullList, char subunit , char domain) {
		int[][] parsingQ3R; 
		TricYeastAlignment alignments = new TricYeastAlignment();
		switch (domain) {
    	case 'B': // Both Equatorial and Middle
    		int[][] newParse1 = {{18,217} , {369,523}};
    		parsingQ3R = newParse1;
    		break;
    	case 'E': // Equatorial 
    		int[][] newParse2 = {{18,149} , {404,523}};
    		parsingQ3R = newParse2;
    		break;
    	case 'M': // Middle 
    		int[][] newParse3 = {{150,217} , {369,403}};
    		parsingQ3R = newParse3;
    		break;
    	case 'A': // Apical 
    		int[][] newParse4 = {{218,368}};
    		parsingQ3R = newParse4;
    		break;
    	default:
    		throw new RuntimeException("Invalid domain letter {E,A}");
		}
		
		int[][] parsing = new int[parsingQ3R.length][2];
		for (int c=0 ; c<parsingQ3R.length ; c++) {
			parsing[c][0] = alignments.getNewResNum('K', parsingQ3R[c][0], subunit);
			//System.out.print(parsing[c][0] + ",");
			parsing[c][1] = alignments.getNewResNum('K', parsingQ3R[c][1], subunit);
			//System.out.print(parsing[c][1] + ",");
		}	
//		if (domain=='E')
//			System.out.println(parsing[0][0]); 
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
