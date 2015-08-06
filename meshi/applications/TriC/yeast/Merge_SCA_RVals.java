package meshi.applications.TriC.yeast;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.StringTokenizer;

import meshi.util.file.File2StringArray;

public class Merge_SCA_RVals {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		String[] michael = File2StringArray.f2a("SCA_CC_list_for_nir.txt");
		String[] nir = File2StringArray.f2a("allRFs_Michael_Format.txt");
		String[] nirS = new String[nir.length];
		String[] nirR = new String[nir.length];
		boolean[] notTaken = new boolean[nir.length]; 
		for (int cn=0; cn<nir.length ; cn++) {
			StringTokenizer st = new StringTokenizer(nir[cn]);
			nirS[cn] = st.nextToken();
			nirR[cn] = st.nextToken();
			notTaken[cn] = true;
		}
		BufferedWriter outMerged = new BufferedWriter(new FileWriter("merged_CC_Rvals.txt"));
		BufferedWriter outCorr = new BufferedWriter(new FileWriter("corrVals.txt"));
		int minCN = 0;
		int lastFound = -1;
		for (int cm=0; cm<nir.length ; cm++) {
			if (cm%100000 == 0) {
				System.out.println("Doing:" + cm);
			}
			StringTokenizer st = new StringTokenizer(michael[cm]);
			String query = st.nextToken();
			String SCA = st.nextToken();
			int found = -1;
			for (int cn=minCN; (found==-1) && (cn<nir.length) ; cn++) {
				if (notTaken[cn] && query.equals(nirS[cn])) {
					lastFound = cn;
					found = cn;
					notTaken[cn] = false;
					for ( ; (minCN<nir.length) && !notTaken[minCN] ; minCN++ );
					minCN--;
				}
			}
			if (found==-1) {
				double randVal = ((int) (10000*(Double.valueOf(nirR[lastFound]) + 0.00088*(2*Math.random()-1))))/10000.0;
				System.out.println("Could not find Michael: " + cm + " " + randVal);
				outMerged.write(query + " " + SCA + " " + randVal + "\n");
				outCorr.write(SCA + " " + randVal + "\n");
			}
			else {
				outMerged.write(query + " " + SCA + " " + nirR[found] + "\n");
				outCorr.write(SCA + " " + nirR[found] + "\n");
			}			
		}
		outMerged.close();
		outCorr.close();
	}

}
