package alignment;

import meshi.util.file.File2StringArray;

import java.util.StringTokenizer;

public class RvalSequence extends Sequence {
	
	public final double gapOpeningIncentive = 20.0;
	
	public RvalSequence(String RvalFileName) {
		int gapCounter=0;
		String[] Rfile = File2StringArray.f2a(RvalFileName);
		int lastResNum = Integer.MIN_VALUE;
		for (String str : Rfile) {
			StringTokenizer tok = new StringTokenizer(str);
			int resNum = Integer.parseInt(tok.nextToken());
			int resType = Integer.parseInt(tok.nextToken());
			double[] Rfits = new double[20];
			for (int c=0 ; c<20 ; c++) {
				Rfits[c] = Double.parseDouble(tok.nextToken());
			}
			if ((resNum-1 != lastResNum) && (lastResNum != Integer.MIN_VALUE)) {
				((RvalPosition) lastElement()).setGapOpeningScore(gapOpeningIncentive);
				gapCounter++;
			}
			RvalPosition rpos = new RvalPosition(resNum,resType,Rfits);
			add(rpos);
			lastResNum = resNum;
		}
		((RvalPosition) lastElement()).setGapOpeningScore(0.0); // To allow extension at the end for Needleman-Wunch
//		System.out.println("The gap opening reward is: " + gapOpeningIncentive);
//		System.out.println("Number of structural gaps: " + gapCounter);
	}
	
}