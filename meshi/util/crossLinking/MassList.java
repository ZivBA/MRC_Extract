package meshi.util.crossLinking;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.Vector;


public class MassList implements ResidueMasses {
	
	public MassList() {}
	
	private int letterLUT(char seqLet) {
	  	switch (seqLet) {
    	case 'A': 
    		return 0;
    	case 'C': 
    		return 1;
    	case 'D': 
    		return 2;
    	case 'E': 
    		return 3;
    	case 'F': 
    		return 4;
    	case 'G': 
    		return 5;
    	case 'H': 
    		return 6;
    	case 'I': 
    		return 7;
    	case 'K': 
    		return 8;
    	case 'L': 
    		return 9;
    	case 'M': 
    		return 10;
    	case 'N': 
    		return 11;
    	case 'P': 
    		return 12;
    	case 'Q': 
    		return 13;
    	case 'R': 
    		return 14;
    	case 'S': 
    		return 15;
    	case 'T': 
    		return 16;
    	case 'V': 
    		return 17;
    	case 'W': 
    		return 18;
    	case 'Y': 
    		return 19;
    	case 'Z': 
    		return 20;
    	default:
    		throw new RuntimeException("Invalid sequene letter");
    	}
	}
	
	
	public double getMass(String seq) {
		double mass = 0.0;
		for (int c=0 ; c<seq.length() ; c++) {
			mass += MW_AA[letterLUT(seq.charAt(c))];
		}
		mass += MW_H2O;
		return mass;				
	}
	
	public double getMass(Digestion digest , int pepInd) {
		return getMass(digest.getSeq(pepInd));
	}

	public void printMassList(Digestion digest) {
		int counter = 0;
		for (int c=0 ; c<digest.getNumOfPeptides() ; c++) {
			double[] masses = getMasses(digest,c);
			for (int m=0 ; m<masses.length ; m++) {
				counter++;
//				System.out.println(counter + ". " + masses[m]);
				System.out.println(masses[m] + " " + (c+1));
			}
		}
	}

	public double[] getMasses(Digestion digest , int pepInd) {
		return getMasses(digest.getSeq(pepInd));
	}

	public double[] getMasses(String seq) {
		Vector<Double> masses = new Vector<Double>();
		masses.add(new Double( MW_H2O ));
//		masses.add(new Double( MW_H2O + 21.981944));
		for (int c=0 ; c<seq.length() ; c++) {
			int currentLengthOfMassesVec = masses.size();
			for (int m=1 ; m<MW_AA_MODIfICATIONS[letterLUT(seq.charAt(c))].length ; m++) {
				for (int massC=0 ; massC<currentLengthOfMassesVec ; massC++) {
					masses.add(new Double( 
							masses.get(massC).doubleValue() + MW_AA_MODIfICATIONS[letterLUT(seq.charAt(c))][m]  
							));
				}
			}
			for (int massC=0 ; massC<currentLengthOfMassesVec ; massC++) {
				masses.set(massC, new Double( 
						masses.get(massC).doubleValue() + MW_AA_MODIfICATIONS[letterLUT(seq.charAt(c))][0]  
						));
			}
		}
		
		double[] outVec = new double[masses.size()];
		for (int c=0 ; c<masses.size() ; c++) {
			outVec[c] = masses.get(c).doubleValue();
		}
		return outVec;				
	}

	public void printXLMassList(Digestion digest1 , Digestion digest2, String outFilename) {
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFilename));
			int counter = 0;
			for (int d1=0 ; d1<digest1.getNumOfPeptides() ; d1++) {
				double[] masses1 = getMasses(digest1,d1);
				for (int d2=0 ; d2<digest2.getNumOfPeptides() ; d2++) {
					double[] masses2 = getMasses(digest2,d2);
					for (int m1=0 ; m1<masses1.length ; m1++) {
						for (int m2=0 ; m2<masses2.length ; m2++) {
							counter++;
							bw.write((masses1[m1]+masses2[m2]+MW_CL) + " " + (d1+1) + " " + (d2+1) + "\n");
							System.out.println((masses1[m1]+masses2[m2]+MW_CL) + " " + (d1+1) + " " + (d2+1));						
							//						System.out.println(counter + ". " + (masses1[m1]+masses2[m2]+MW_CL) + " " + (d1+1) + " " + (d2+1));						
						}					
					}
				}
			}
			bw.close();
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
    	}
	}

	
}
