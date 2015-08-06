package meshi.applications.TriC.yeast;

import meshi.util.file.File2StringArray;

public class MakingTricAlignmentTableS3 {

	public static void main(String[] args) {
		String[] clustal = File2StringArray.f2a("C:\\Users\\Nir\\TRiC\\Organisms\\Alignment_for_Fig_S3");
		String[] headers = {"1Q3R=A  ","YeastA  ","YeastB  ","YeastG  ","YeastD  ","YeastE  ","YeastZ  ","YeastH  ","YeastQ  "};
		String[] names = {"1Q3R   ","YeastA ","YeastB ","YeastG ","YeastD ","YeastE ","YeastZ ","YeastH ","YeastQ "};
		
		int[] counters = new int[headers.length];
		for (int c=0 ;  c<10 ; c++) {
			int firstInd = 8 + c*65;
			int lastInd = 8 + (c+1)*65;
			if (c==9) {
				lastInd = 651;
			}
			for (int seq=0; seq<headers.length ; seq++) {
				int counter=0;
				for (int tmp=0 ; tmp<clustal[seq].substring(firstInd, lastInd).length() ; tmp++) {
					if (clustal[seq].substring(firstInd, lastInd).charAt(tmp)!='-')
						counter++;
				}
				counters[seq] += counter;
				System.out.println(names[seq] + clustal[seq].substring(firstInd, lastInd) + "  " +counters[seq]);				
			}
			String core = "";
			String concensus = "";
			for (int tmp=0 ; tmp<clustal[0].substring(firstInd, lastInd).length() ; tmp++) {
				boolean same = true;
				boolean coreTrue = true;
				for (int seq=0; seq<headers.length ; seq++) {
					if (clustal[0].substring(firstInd, lastInd).charAt(tmp) != clustal[seq].substring(firstInd, lastInd).charAt(tmp)) {
						same = false;
					}
					if (clustal[seq].substring(firstInd, lastInd).charAt(tmp) == '-') {
						coreTrue = false;
					}
				}
				if (same) 
					concensus += "*";
				else
					concensus += " ";
				if (coreTrue) 
					core += "^";
				else
					core += " ";
			}
			System.out.println("       " + core);
			System.out.println("       " + concensus + "\n");
		}

	}

}
