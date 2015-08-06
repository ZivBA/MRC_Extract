package programs;

import java.util.Arrays;
import java.util.StringTokenizer;

import meshi.util.file.File2StringArray;

/**
 * Running:
 * 
 * java EnrichmentPhase1ParamaterSearch <list of data files> 
 * 
 * @author Nir
 *
 */

public class EnrichmentPhase1ParamaterSearch {
	
	public static void main(String[] args) {
		System.out.println("Reading files from: " + args[0].trim());

		// Enrichment percentage
		double enrichPercent = 0.1;
		
		// Columns - Parameters matches
		int column1 = 5; // Centroid
		int column2 = 6; // HB SR
		int column3 = 7; //HB LR 
		int column4 = 8;  // EV
		int column5 = 9;  // PROP
		int column6 = 10;  // Clust
		int columnRMS = 4;
		
		// parameter ranges
		double[] range1 = {-0.2, 0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0};
		double[] range2 = {-0.2, 0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0};
		double[] range3 = {1.0};
		double[] range4 = {-0.2, 0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0};
		double[] range5 = {-0.2, 0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0};
		double[] range6 = {0.2, 0, -0.05, -0.1, -0.25, -0.5, -1.0, -2.0, -4.0};
		
		
		// Reading data
		String[] listOfDataFiles = File2StringArray.f2a(args[0].trim());
		double[][][] data = new double[listOfDataFiles.length][][];
		double[] pertileRMS = new double[data.length];
		double[] pertile1 = new double[data.length];
		double[] pertile2 = new double[data.length];
		double[] pertile3 = new double[data.length];
		double[] pertile4 = new double[data.length];
		double[] pertile5 = new double[data.length];
		double[] pertile6 = new double[data.length];
		String token;
		for (int loop=0; loop<data.length ; loop++) {
			String[] tmpData = File2StringArray.f2a(listOfDataFiles[loop]);
			data[loop] = new double[9][tmpData.length];
			for (int decoy=0 ; decoy<tmpData.length ; decoy++) {
				StringTokenizer st = new StringTokenizer(tmpData[decoy]);
				for (int column=1 ; st.hasMoreTokens() ; column++) {
					token = st.nextToken();
					if (column==columnRMS)
						data[loop][0][decoy] = (new Double(token)).doubleValue();
					if (column==column1)
						data[loop][1][decoy] = (new Double(token)).doubleValue();
					if (column==column2)
						data[loop][2][decoy] = (new Double(token)).doubleValue();
					if (column==column3)
						data[loop][3][decoy] = (new Double(token)).doubleValue();
					if (column==column4)
						data[loop][4][decoy] = (new Double(token)).doubleValue();
					if (column==column5)
						data[loop][5][decoy] = (new Double(token)).doubleValue();
					if (column==column6)
						data[loop][6][decoy] = (new Double(token)).doubleValue();
				}
			}
			pertileRMS[loop] = getPercentile(data[loop][0],enrichPercent);
			pertile1[loop] = getPercentile(data[loop][1],enrichPercent);
			pertile2[loop] = getPercentile(data[loop][2],enrichPercent);
			pertile3[loop] = getPercentile(data[loop][3],enrichPercent);
			pertile4[loop] = getPercentile(data[loop][4],enrichPercent);
			pertile5[loop] = getPercentile(data[loop][5],enrichPercent);
			pertile6[loop] = getPercentile(data[loop][6],enrichPercent);
		}
		
		// Do the parameter search
		double bestEnrich = -9999.0;
		double totalEnrich = 0;
		double loopEnrich = 0;
		double energyPercentile = 0;
		double best1 = -9999;
		double best2 = -9999;
		double best3 = -9999;
		double best4 = -9999;
		double best5 = -9999;
		double best6 = -9999;
		for (int c3=0 ; c3<range3.length ; c3++)
		for (int c1=0 ; c1<range1.length ; c1++)
		for (int c2=0 ; c2<range2.length ; c2++)
		for (int c4=0 ; c4<range4.length ; c4++)
		for (int c5=0 ; c5<range5.length ; c5++) 
		for (int c6=0 ; c6<range6.length ; c6++) {
			// Calculating the enrichment
			totalEnrich = 0.0;
			for (int loop=0 ; loop<data.length ; loop++) {
				for (int decoy=0 ; decoy<data[loop][0].length ; decoy++) {
					data[loop][7][decoy] = range1[c1]*data[loop][1][decoy] + 
					range2[c2]*data[loop][2][decoy] + 
					range3[c3]*data[loop][3][decoy] + 
					range4[c4]*data[loop][4][decoy] + 
					range5[c5]*data[loop][5][decoy] + 
					range6[c6]*data[loop][6][decoy];
					data[loop][8][decoy] = data[loop][7][decoy];
				}
				Arrays.sort(data[loop][8]);
				energyPercentile = data[loop][8][(int) (enrichPercent*data[loop][8].length)];
				loopEnrich = 0;
				for (int decoy=0 ; decoy<data[loop][0].length ; decoy++) { 
					if ((data[loop][7][decoy]<energyPercentile) && (data[loop][0][decoy]<pertileRMS[loop]))
						loopEnrich++;
				}
				loopEnrich = loopEnrich/(enrichPercent*enrichPercent*data[loop][0].length);
				totalEnrich += loopEnrich;
			}
			totalEnrich = totalEnrich/data.length;
			// Is this better?
			if (totalEnrich>bestEnrich) {
				best1 = range1[c1];
				best2 = range2[c2];
				best3 = range3[c3];
				best4 = range4[c4];
				best5 = range5[c5];
				best6 = range6[c6];
				bestEnrich = totalEnrich;
				System.out.println("Best: " + bestEnrich + " " + best1 + " " + best2 + " " + best3 + " " + best4 + " " + best5 + " " + best6);
			}
		}		
	}
	
	public static double getPercentile(double[] ar, double pertile) {
		double[] copy_ar = new double[ar.length];
		for (int c=0 ; c<ar.length ; c++)
			copy_ar[c] = ar[c];
		Arrays.sort(copy_ar);
		return copy_ar[(int) (pertile*copy_ar.length)];
	}

}
