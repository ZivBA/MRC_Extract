package meshi.util.Isite;

import java.util.StringTokenizer;
import java.util.Vector;

public class Imotif {
	
	private int clusterID = -1;
	private int motifLength = -1; // The profile is always 2*overhang positions more
	private int overhang = -1;
	private double mdaCut = -1;
	private double dmeCut = -1;
	private double pseudoCount = -1;
	private double confidenceFitParam1 = -1;
	private double confidenceFitParam2 = -1;
	private double zScoreParam1 = -1;
	private double zScoreParam2 = -1;
//	private double confidenceFitParam3 = -1;
//	private double confidenceFitParam4 = -1;
	private double Wcovar = -1;
	private String paradigm = "";
	private double[][] angles = null;
	private double[][] profile = null;
//	private double[][][][] covar = null;
	                 
	public Imotif(Vector<String> libString) {
		int lineCounter = 0;
		StringTokenizer st = new StringTokenizer(libString.get(0)); // Parsing cluster line
		st.nextToken();
		clusterID = Integer.parseInt(st.nextToken());
		motifLength = Integer.parseInt(st.nextToken());
		st.nextToken();
		overhang = Integer.parseInt(st.nextToken());
		st = new StringTokenizer(libString.get(3)); // Parsing MDA  
		st.nextToken();
		mdaCut = Math.PI/180.0*Double.parseDouble(st.nextToken());
		st = new StringTokenizer(libString.get(4)); // Parsing DME  
		st.nextToken();
		dmeCut = Double.parseDouble(st.nextToken());
		st = new StringTokenizer(libString.get(5)); // Parsing Pseudo Counts  
		st.nextToken();
		pseudoCount = Double.parseDouble(st.nextToken());
//		st = new StringTokenizer(libString.get(6)); // Parsing confidence params  
//		st.nextToken();
//		confidenceFitParam1 = Double.parseDouble(st.nextToken());
//		confidenceFitParam2 = Double.parseDouble(st.nextToken());
//		confidenceFitParam3 = Double.parseDouble(st.nextToken());
//		confidenceFitParam4 = Double.parseDouble(st.nextToken());
		st = new StringTokenizer(libString.get(7)); // Parsing Wcovar  
		st.nextToken();
		Wcovar = Double.parseDouble(st.nextToken());
		st = new StringTokenizer(libString.get(8)); // Parsing paradigme  
		st.nextToken();
		paradigm = st.nextToken() + " " + st.nextToken() + " " + st.nextToken();
		angles = new double[motifLength][3];
		for (lineCounter=10 ; lineCounter<(10+motifLength) ; lineCounter++) {  // Parsing the ANGLES
			st = new StringTokenizer(libString.get(lineCounter));
			st.nextToken();
			angles[lineCounter-10][0] = Math.PI/180.0*Double.parseDouble(st.nextToken());
			angles[lineCounter-10][1] = Math.PI/180.0*Double.parseDouble(st.nextToken());
			angles[lineCounter-10][2] = Math.PI/180.0*Double.parseDouble(st.nextToken());
		}
		lineCounter++;
		profile = new double[motifLength+2*overhang][20];
		for (int c=0; c <(motifLength+2*overhang) ; c++) {  // Parsing the profile
			st = new StringTokenizer(libString.get(lineCounter+c));
			st.nextToken();
			for (int aa=0 ; aa<20 ; aa++) {
				profile[c][aa] = Double.parseDouble(st.nextToken());
			}			
		}
//		lineCounter += (motifLength+2*overhang+1);
//		covar = new double[motifLength][motifLength][9][9];
//		for (int c=0; c<motifLength ; c++) {  // Parsing the covar
//			for (int d=c+1; d<motifLength ; d++) {
//				lineCounter++;
//				for (int t1=0 ; t1<9 ; t1++) {
//					st = new StringTokenizer(libString.get(lineCounter));
//					for (int t2=0 ; t2<9 ; t2++) {
//						covar[c][d][t1][t2] = Double.parseDouble(st.nextToken());						
//					}
//					lineCounter++;
//				}
//			}
//		}		
	}


	public int getClusterID() {
		return clusterID;
	}

	public int getMotifLength() {
		return motifLength;
	}
	
	public int getOverhang() {
		return overhang;
	}

	public String getParadigm() {
		return paradigm;
	}
	
	public double getPhi(int ind) {
		return angles[ind][0];
	}
	public double getPsi(int ind) {
		return angles[ind][1];
	}
	public double getOmg(int ind) {
		return angles[ind][2];
	}
	
	public double getProf(int pos, int aa) {
		return profile[pos][aa];
	}

	public double getDMAcut() {
		return mdaCut;
	}

	public double getDMEcut() {
		return dmeCut;
	}
	
	public double getWcovar() {
		return Wcovar;
	}
	
	public void setConfidenceParams(double a, double b) {
		confidenceFitParam1 = a;
		confidenceFitParam2 = b;
	}

	public double confidenceA() {
		return confidenceFitParam1;
	}

	public double confidenceB() {
		return confidenceFitParam2;
	}

	public void setZscoreParams(double a, double b) {
		zScoreParam1 = a;
		zScoreParam2 = b;
	}

	public double zScoreMean() {
		return zScoreParam1;
	}

	public double zScoreStd() {
		return zScoreParam2;
	}

}
