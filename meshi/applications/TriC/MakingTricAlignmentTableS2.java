package meshi.applications.TriC;

import meshi.util.file.File2StringArray;

public class MakingTricAlignmentTableS2 {

	public static void main(String[] args) {
		String[] clustal = File2StringArray.f2a("C:\\Users\\Nir\\TRiC\\Organisms\\Alignment_for_Fig_S2");
		String[] headers = {"1Q3R=A  ","1A6D=A  ","1A6D=B  ","mmCPNA  ","TRiC=A  ","TRiC=B  ","TRiC=G  ","TRiC=D  ","TRiC=E  ","TRiC=Z  ","TRiC=W  ","TRiC=H  ","TRiC=Q  "};
		String[] names = {"1Q3R    ","1A6D_A  ","1A6D_B  ","mmCPN   ","BovA    ","BovB    ","BovG    ","BovD    ","BovE    ","BovZ    ","BovW    ","BovH    ","BovQ    "};
		
		int[] counters = new int[headers.length];
		for (int c=0 ;  c<10 ; c++) {
			int firstInd = 8 + c*65;
			int lastInd = 8 + (c+1)*65;
			if (c==9) {
				lastInd = 613;
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
			String concensus = "";
			for (int tmp=0 ; tmp<clustal[0].substring(firstInd, lastInd).length() ; tmp++) {
				boolean same = true;
				for (int seq=0; seq<headers.length ; seq++) {
					if (clustal[0].substring(firstInd, lastInd).charAt(tmp) != clustal[seq].substring(firstInd, lastInd).charAt(tmp)) {
						same = false;
					}
				}
				if (same) 
					concensus += "*";
				else
					concensus += " ";
			}
			System.out.println("        " + concensus + "\n");
		}

	}

}
