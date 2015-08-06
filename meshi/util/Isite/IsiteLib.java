package meshi.util.Isite;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.StringTokenizer;
import java.util.Vector;

import meshi.applications.corpus.Corpus;
import meshi.applications.loopBuilding.AbstractLoopBuilder;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.Residues;
import meshi.util.rotamericTools.RotamericTools;

public class IsiteLib extends Vector<Imotif> implements Residues {
	
	private final int motifDepthInCorpus = 5;
	private Corpus corpus = null;
	private int[][] motifInCorpus = null;
	private int[][] motifPosInCorpus = null;
	private double[][] motifScoresInCorpus = null;
	private int[] loadedSeq = null;
	private int[][] posInSeq = null;
	private double[][] scoresInSeq = null;
	private double[][] phiPsiInSeq = null;
	private int[] sortedIndicesByDecendingMotifLength = null;
	public final static double veryHighEnergy = 100000000000.0;

	//	private final double deltaConfToOverride = 0.25;
	
	public IsiteLib(String fileName) {
		System.out.println("Reading I-site lib: " + fileName);
		try{
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			String line = br.readLine();
			while (!line.startsWith("CLUSTER ")) {
				line = br.readLine();
			}
			do {
				//System.out.println("Adding new library:" + line);
				Vector<String> libString = new Vector<String>();
				libString.add(line);
				line = br.readLine();
				while ((line!=null) && (!line.startsWith("CLUSTER "))) {
					libString.add(line);
					line = br.readLine();
				}
				add(new Imotif(libString));
				if (lastElement().getClusterID()!=confCurvIsite51[size()-1][0]) {
					throw new RuntimeException("A mismatch with the confidence paremeters calculated by Nir");
				}
				else {
					lastElement().setConfidenceParams(confCurvIsite51[size()-1][1], confCurvIsite51[size()-1][2]);
				}
				if (lastElement().getClusterID()!=zScoresIsite51[size()-1][0]) {
					throw new RuntimeException("A mismatch with the Z-scores paremeters calculated by Nir");
				}
				else {
					lastElement().setZscoreParams(zScoresIsite51[size()-1][1], zScoresIsite51[size()-1][2]);
				}
				//System.out.println((size()-1) + " " + lastElement().getClusterID());
			} while (line!=null);
		    br.close();
		}
		catch(Exception e) {
		    throw new RuntimeException(e.getMessage());
		}
		double[] motifLengths = new double[size()];
		for (int c=0 ; c<size() ; c++) {
			motifLengths[c] = -get(c).getMotifLength();
		}
		sortedIndicesByDecendingMotifLength = AbstractLoopBuilder.findTopMinArray(motifLengths, size(), Double.MAX_VALUE);
	}
	
	
	/**
	 * Assign a corpus to the library
	 * @param corpus
	 */
	public void setCorpus(Corpus corpus) {
		this.corpus = corpus;
	}
	
	public int findCenterOfMotif(Imotif motif) {
		int overhang = motif.getOverhang();
		int motifLength = motif.getMotifLength();
		corpus.buildUngappedArray(motifLength + 2*overhang);
		double minMDA = Double.MAX_VALUE;
		int bestFrag = -1;
		for (int c=0 ; c<corpus.ungapped.length ; c++) {
			double mda = calcMDA(motif,corpus,corpus.ungapped[c]+overhang);
			if (mda<minMDA) {
				minMDA = mda;
				bestFrag = corpus.ungapped[c]+overhang;
//				System.out.println(minMDA + " " + motif.getClusterID() + " " + motif.getParadigm() + 
//						"            " + corpus.proteinNames[corpus.protInd[bestFrag]] + " " + 
//						corpus.resNum[bestFrag] + "       score: " + 
//						calcScore(motif, corpus, bestFrag-overhang));		
			}
		}
//		System.out.println(minMDA + " " + motif.getClusterID() + " " + motif.getParadigm() + "            " + corpus.proteinNames[corpus.protInd[bestFrag]] + " " + corpus.resNum[bestFrag] + "   " + motif.getDMAcut());		

		return bestFrag;
	}

	
	public void calcSuccess(int motifInd) {
		Imotif motif = get(motifInd);
		int overhang = motif.getOverhang();
		int motifLength = motif.getMotifLength();
		int center = findCenterOfMotif(motif);
		corpus.buildUngappedArray(motifLength + 2*overhang);
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("success_"+motif.getClusterID()+".dat")));	
			pw.println(motif.getWcovar() + "  " + motif.getMotifLength());
			System.out.println(motif.getWcovar() + "  " + motif.getMotifLength());
			for (int c=0 ; c<corpus.ungapped.length ; c++) {
				if (isFragInMotif(motif, corpus, center, corpus.ungapped[c]+overhang)) {
//					System.out.println(1 + " " + calcScore(motif, corpus, corpus.ungapped[c]));
					pw.println(1 + " " + calcScore(motif, corpus, corpus.ungapped[c]));
				}
				else {
//					System.out.println(0 + " " + calcScore(motif, corpus, corpus.ungapped[c]));
					pw.println(0 + " " + calcScore(motif, corpus, corpus.ungapped[c]));
				}
			}
			pw.close();
		}
		catch (Exception e) {
			throw new RuntimeException(e);
		}

	}

	public static boolean isFragInMotif(Imotif motif, Corpus corpus, int centerInCorpus, int motifStart) {
		double dmaCut = motif.getDMAcut();
		double dmeCut = motif.getDMEcut();
		double mda = calcMDA(motif,corpus,motifStart);
		double dme = calcDME(motif, corpus, centerInCorpus, motifStart);
//		if ((mda<dmaCut) && (dme<dmeCut)) {
		if ((mda<Math.min(Math.PI/6.0+motif.getMotifLength()/15.0*Math.PI/6.0,dmaCut)) && (dme<Math.min(dmeCut,1.5+0.25*(motif.getMotifLength()-8)))) {
			return true;
		}
		else {
			return false;
		}		
	}
	

	public void calcMotifsInCorpus() {
		motifInCorpus = new int[corpus.Nres][motifDepthInCorpus];
		motifPosInCorpus = new int[corpus.Nres][motifDepthInCorpus];
		motifScoresInCorpus = new double[corpus.Nres][motifDepthInCorpus];
		double[] mdas = new double[corpus.Nres];
		double[] dmes = new double[corpus.Nres];
		DecimalFormat fmt = new DecimalFormat("0.##"); 
		for (int c=0 ; c<corpus.Nres ; c++) {
			for (int depth=0 ; depth<motifDepthInCorpus ; depth++) {
				motifInCorpus[c][depth] = -1;
				motifPosInCorpus[c][depth] = -1;
				motifScoresInCorpus[c][depth] = -10000;
			}
		}
		for (int counter=0 ; counter<size() ; counter++) {
			int motifInd = sortedIndicesByDecendingMotifLength[counter];
			Imotif motif = get(motifInd);
			System.out.print("Calculating: " + motif.getClusterID());
			int overhang = motif.getOverhang();
			int motifLength = motif.getMotifLength();
			int center = findCenterOfMotif(motif);
			corpus.buildUngappedArray(motifLength + 2*overhang);
			int overturnNum = 0;
			int assigned = 0;
			for (int c=0 ; c<corpus.ungapped.length ; c++) {
				if (isFragInMotif(motif, corpus, center, corpus.ungapped[c]+overhang)) {
					boolean overturned = false;
//					int overturnLength = -1;
//					int overturnMotif = -1;
//					double overturnScore = -1;
//					double overturnmdas = -1;
//					double overturndmes = -1;
					double score = calcScore(motif, corpus, corpus.ungapped[c]);
					for (int pos=0 ; pos<(motifLength + 2*overhang) ; pos++) {
						int overturnInd = findWhereToOverTurn(motifScoresInCorpus[corpus.ungapped[c]+pos]);
						if (score>motifScoresInCorpus[corpus.ungapped[c]+pos][overturnInd]) {
//						if ((score>motifScoresInCorpus[corpus.ungapped[c]+pos][overturnInd]+deltaConfToOverride) || 
//								((score>motifScoresInCorpus[corpus.ungapped[c]+pos][overturnInd]) && 
//										(motifInCorpus[corpus.ungapped[c]+pos][overturnInd]==motifInd))) {
							if ((motifScoresInCorpus[corpus.ungapped[c]+pos][overturnInd]>-10000) &&
									((motif.getMotifLength()+1)<get(motifInCorpus[corpus.ungapped[c]+pos][overturnInd]).getMotifLength())) {
								overturned = true;
//								overturnLength = get(motifInCorpus[corpus.ungapped[c]+pos]).getMotifLength();
//								overturnMotif = get(motifInCorpus[corpus.ungapped[c]+pos]).getClusterID();
//								overturnScore = motifScoresInCorpus[corpus.ungapped[c]+pos];
//								overturnmdas = mdas[corpus.ungapped[c]+pos];
//								overturndmes = dmes[corpus.ungapped[c]+pos];
								overturnNum++;
							}
							motifScoresInCorpus[corpus.ungapped[c]+pos][overturnInd] = score;
							mdas[corpus.ungapped[c]+pos] = calcMDA(motif,corpus,corpus.ungapped[c]+overhang);
							dmes[corpus.ungapped[c]+pos] = calcDME(motif,corpus,center,corpus.ungapped[c]+overhang);
							motifPosInCorpus[corpus.ungapped[c]+pos][overturnInd] = pos;
							motifInCorpus[corpus.ungapped[c]+pos][overturnInd] = motifInd;
							assigned++;
						}
					}
					if (overturned) {
//						System.out.println("The " + motif.getClusterID() + " of length " + motif.getMotifLength() +
//								" [" + fmt.format(score) + " " + fmt.format(calcMDA(motif,corpus,corpus.ungapped[c]+overhang)) +
//								" " + calcDME(motif,corpus,center,corpus.ungapped[c]+overhang) + "]   " + 
//								" overturned " + overturnMotif + " of length " + overturnLength + 
//								" [" + fmt.format(overturnScore) + " " + fmt.format(overturnmdas) + " " + fmt.format(overturndmes) + "]");
					}
				}
			}
			System.out.println(" assigned: " + assigned + " " + assigned*100.0/corpus.Nres + "%     overturn: " + overturnNum);
		}
		
		// Calcultaing statistics:
		int coverage = 0;
		for (int c=0 ; c<corpus.Nres ; c++) {
			if (motifInCorpus[c][0]!=-1) {
				coverage++;
			}
		}
		System.out.println("The coverage is: " + coverage*100.0/corpus.Nres);
		
		
		/// Saving it
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("motifAssignment.dat")));	
			for (int c=0 ; c<corpus.Nres ; c++) {
				for (int depth=0 ; depth<motifDepthInCorpus ; depth++) {
					if (motifInCorpus[c][depth]==-1) {
						pw.print(motifInCorpus[c][depth] + " " + motifPosInCorpus[c][depth] + " " + fmt.format(motifScoresInCorpus[c][depth]) + " 00000   ");
					}
					else {
						pw.print(motifInCorpus[c][depth] + " " + motifPosInCorpus[c][depth] + " " + fmt.format(motifScoresInCorpus[c][depth]) + " " + get(motifInCorpus[c][depth]).getClusterID() + "   ");
					}
				}
				pw.println();
			}
			pw.close();
		}
		catch (Exception e) {
			throw new RuntimeException(e);
		}	
	}
	
	private static int findWhereToOverTurn(double[] scores) {
		double min = Double.MAX_VALUE;
		int whereMin = -1;
		for (int answer=0; answer<scores.length ; answer++) {
			if (scores[answer]==-10000) {
				return answer;
			}
			if (scores[answer]<min) {
				whereMin = answer;
				min = scores[answer];
			}
		}
		return whereMin;
	}
	
	
	public void loadMotifInCorpus(String motifFileName) {
		motifInCorpus = new int[corpus.Nres][motifDepthInCorpus];
		motifPosInCorpus = new int[corpus.Nres][motifDepthInCorpus];
		motifScoresInCorpus = new double[corpus.Nres][motifDepthInCorpus];
		try{
			BufferedReader br = new BufferedReader(new FileReader(motifFileName));
			String line;
			for (int c=0 ; c<corpus.Nres ; c++) {
				line = br.readLine();
				StringTokenizer st = new StringTokenizer(line);
				for (int depth=0 ; depth<motifDepthInCorpus ; depth++) {
					motifInCorpus[c][depth] = Integer.parseInt(st.nextToken());
					motifPosInCorpus[c][depth] = Integer.parseInt(st.nextToken());				
					motifScoresInCorpus[c][depth] = Double.parseDouble(st.nextToken());
					st.nextToken(); // The clusterID
				}
			}
		    br.close();
		}
		catch(Exception e) {
		    throw new RuntimeException(e.getMessage());
		}		
	}
	
	public void loadAndAnalyzeSequenceIntoLib(int[] seq) {
		// Loading sequence
		loadedSeq = new int[seq.length];
		for (int c=0; c<seq.length ; c++) {
			if ((seq[c]<0) || (seq[c]>19)) { // Padding with alanines.
				loadedSeq[c] = 0;				
			}
			else {
				loadedSeq[c] = seq[c];
			}
		}
		// Analyzing
		posInSeq = new int[size()][loadedSeq.length];
		scoresInSeq = new double[size()][loadedSeq.length];
		for (int motifInd=0 ; motifInd<size() ; motifInd++) {
			Imotif motif = get(motifInd);
			//System.out.println("Calculating: " + motif.getClusterID());
			int overhang = motif.getOverhang();
			int motifLength = motif.getMotifLength();
			int[] tmpSeq = new int[2*overhang+motifLength];
			for (int c=0; c<loadedSeq.length ; c++) { 
				posInSeq[motifInd][c] = -1;
				scoresInSeq[motifInd][c] = -10000;
			}
			for (int c=0 ; c<loadedSeq.length-2*overhang-motifLength+1 ; c++) {
				for (int pos=0 ; pos<2*overhang+motifLength ; pos++) {
					tmpSeq[pos] = loadedSeq[c+pos];
				}
				double score = calcScoreFromSeq(motif, tmpSeq);
				for (int pos=0 ; pos<(motifLength + 2*overhang) ; pos++) {
					if (score>scoresInSeq[motifInd][c+pos]) {
						posInSeq[motifInd][c+pos] = pos;
						scoresInSeq[motifInd][c+pos] = score;
					}
				}
			}
		}		
	}
	
	
	public void loadAndAnalyzeBackboneIntoLib(Protein prot) {
		// Calculating phi,psi
		double [][] pp = RotamericTools.phipsi(prot, new DistanceMatrix(prot.atoms(), 1.0, 0.1, 4));
		phiPsiInSeq = new double [loadedSeq.length][2]; 
		for (int c=0 ;c<loadedSeq.length ; c++) {
			if (pp[c]!=null) {
				phiPsiInSeq[c][0] = pp[c][0];
				phiPsiInSeq[c][1] = pp[c][1];
			}
			else {
				phiPsiInSeq[c][0] = 0.0;
				phiPsiInSeq[c][1] = 0.0;				
			}			
		}
	}

	
	public double[] negativeScoresOfMotifInSeqPos(int pos) {
		double[] answer = new double[size()];
		for (int motifInd=0 ; motifInd<size() ; motifInd++) {
			answer[motifInd] = -scoresInSeq[motifInd][pos];
		}
		return answer;
	}
	
	public int bestMotifInSeqPos(int pos) {
		DecimalFormat fmt = new DecimalFormat("0.##"); 
		double[] scores = negativeScoresOfMotifInSeqPos(pos);
		int[] sortedIndices = AbstractLoopBuilder.findTopMinArray(scores, scores.length, Double.MAX_VALUE);
		System.out.println(pos + " " + Residue.nameOneLetter(loadedSeq[pos]) + " " + 
				get(sortedIndices[0]).getClusterID() + " " + 
				posInSeq[sortedIndices[0]][pos] + " " + fmt.format(scoresInSeq[sortedIndices[0]][pos]) + " " +
				fmt.format(calcConfFromScore(get(sortedIndices[0]), scoresInSeq[sortedIndices[0]][pos])));
		return sortedIndices[0];
	}
	
	
	public static double calcScore(Imotif motif,Corpus corpus, int seqStart) {
		int[] seq = new int[motif.getMotifLength() + 2*motif.getOverhang()];
		for (int pos=0 ; pos<(motif.getMotifLength() + 2*motif.getOverhang()) ; pos++) {
			seq[pos] = corpus.resType[seqStart+pos];
		}
		return calcScoreFromSeq(motif, seq);		
	}

// The complex form:
// -----------------	
//	public static double calcScoreFromSeq(Imotif motif, int[] seq) {
//		double pScore = -1;
//		double[] bkProb = {0.0749, 0.0182, 0.0522, 0.0626, 0.0391, 0.071, 0.0223, 0.0545, 0.0582, 0.0906, 0.0227, 0.0453, 0.0512, 0.0411, 0.0522, 0.0734, 0.0596, 0.0648, 0.0132, 0.0325};
////		double[] bkProbOn1 = new double[20];
////		double[] bkProbOn0 = new double[20];
//		double alpha = 0.5;
//		double mult = 4.0;
////		for (int aa=0 ; aa<20 ; aa++) {
////			bkProbOn1[aa] = Math.log((2*bkProb[aa]+alpha*bkProb[aa])/((1+alpha)*bkProb[aa]));
////			bkProbOn0[aa] = Math.log((bkProb[aa]+alpha*bkProb[aa])/((1+alpha)*bkProb[aa]));
////		}
//		int overhang = motif.getOverhang();
//		int motifLength = motif.getMotifLength();
//		for (int pos=0 ; pos<(motifLength+2*overhang) ; pos++) {
//			for (int aa=0 ; aa<20 ; aa++) {
//				if (seq[pos]==aa) {
//					pScore += Math.log((mult*bkProb[aa]+alpha*bkProb[aa])/((1+alpha)*bkProb[aa]))
//					       *Math.log((motif.getProf(pos,aa)+alpha*bkProb[aa])/((1+alpha)*bkProb[aa]));
//				}
//				else {
//					pScore += Math.log((bkProb[aa]-bkProb[aa]*mult*bkProb[seq[pos]]+alpha*bkProb[aa])/((1+alpha)*bkProb[aa]))
//					    *Math.log((motif.getProf(pos,aa)+alpha*bkProb[aa])/((1+alpha)*bkProb[aa]));					
//				}
//			}
//		}
//		return pScore;	
//	}

	
	public static double calcScoreFromSeq(Imotif motif, int[] seq) {
		double pScore = 0.0;
		double[] bkProb = {0.0749, 0.0182, 0.0522, 0.0626, 0.0391, 0.071, 0.0223, 0.0545, 0.0582, 0.0906, 0.0227, 0.0453, 0.0512, 0.0411, 0.0522, 0.0734, 0.0596, 0.0648, 0.0132, 0.0325};
		double alpha = 0.15;
		int overhang = motif.getOverhang();
		int motifLength = motif.getMotifLength();
		for (int pos=0 ; pos<(motifLength+2*overhang) ; pos++) {
			pScore += Math.log((motif.getProf(pos,seq[pos])+alpha*bkProb[seq[pos]])/((1+alpha)*bkProb[seq[pos]]));
		}		
//		double zScore = (pScore-motif.zScoreMean())/motif.zScoreStd();
//		return zScore;
//		double confScore = 1/(1+Math.exp(motif.confidenceA()*pScore+motif.confidenceB())) + 0.0*motif.getMotifLength();
//		return confScore;
		return pScore;
	}

	public static double calcConfFromScore(Imotif motif, double pScore) {
		double confScore = 1/(1+Math.exp(motif.confidenceA()*pScore+motif.confidenceB())) + 0.0*motif.getMotifLength();
		return confScore;
	}

	
	public static double calcMDA(Imotif motif,Corpus corpus, int anglesStart) {
		double maxDMA = -1;
		for (int c=0 ; c<motif.getMotifLength() ; c++) {
			if (torsionDiff(motif.getPhi(c), corpus.torsions[anglesStart+c][1])>maxDMA) {
				maxDMA = torsionDiff(motif.getPhi(c), corpus.torsions[anglesStart+c][1]);
			}
			if (torsionDiff(motif.getPsi(c), corpus.torsions[anglesStart+c][2])>maxDMA) {
				maxDMA = torsionDiff(motif.getPsi(c), corpus.torsions[anglesStart+c][2]);
			}
			if (torsionDiff(motif.getOmg(c), corpus.torsions[anglesStart+c+1][0])>maxDMA) {
				maxDMA = torsionDiff(motif.getOmg(c), corpus.torsions[anglesStart+c+1][0]);
			}				
		}
		return maxDMA;
	}

	
	public static double calcDME(Imotif motif, Corpus corpus, int center, int fragStart) {
		return corpus.calcRmsBetweenStruct(center, fragStart, motif.getMotifLength(), 1, motif.getMotifLength());
	}

	public static double torsionDiff(double tor1, double tor2) {
		if (Math.abs(tor2-tor1)>Math.PI) {
			return 2*Math.PI - Math.abs(tor2-tor1);
		}
		return Math.abs(tor2-tor1);
	}
	
	public int findSeqStartInProt(int[] seq) {
		for (int c=0 ; c<loadedSeq.length-seq.length+1 ; c++) {
			boolean perfectMatch = true;
			for (int cc=0 ; cc<seq.length ; cc++) {
				if (loadedSeq[c+cc]!=seq[cc]) {
					perfectMatch=false;
				}
			}
			if (perfectMatch) {
				return c;
			}
		}
		return -9999999;
	}
	
	public void setThreadingEnergies(int[] seq, int fragStartInSeq, int fragEndInSeq, double[] energies) {
		int tolerance = 0; // How many residues in the sequence need NOT be covered by the (true==angles) motif. 
		int seqStartInProt = findSeqStartInProt(seq);
		System.out.println("Full sequence starts is residue: " + seqStartInProt);
		System.out.println("True fragment starts is residue: " + (seqStartInProt+fragStartInSeq));
		for (int c=0 ; c<corpus.ungapped.length ; c++) {  // reseting the energies
			energies[c] = veryHighEnergy;
		}
		for (int c=0 ; c<corpus.ungapped.length ; c++) {  // Now actually finding fits
			for (int depth=0; depth<motifDepthInCorpus ; depth++) {
				if (motifPosInCorpus[corpus.ungapped[c]][depth]!=-1) {
					int motifInd = motifInCorpus[corpus.ungapped[c]][depth];
					int motifPos = motifPosInCorpus[corpus.ungapped[c]][depth];
					int motifLength = get(motifInd).getMotifLength();
					int overhang = get(motifInd).getOverhang();
					double motifScore = motifScoresInCorpus[corpus.ungapped[c]][depth];
					for (int pos=seqStartInProt+fragStartInSeq ; pos<=(seqStartInProt+fragEndInSeq) ; pos++) {
						int queryMotifPos = posInSeq[motifInd][pos];
						double queryScore = scoresInSeq[motifInd][pos];
						int fragLength = fragEndInSeq-fragStartInSeq+1;
						int posInFrag = pos - (seqStartInProt+fragStartInSeq);
						// Next checking if the motif overlap with the fragment within tolerance
						if (((queryMotifPos-overhang-posInFrag)>=-(tolerance)) && (((motifLength-(queryMotifPos-overhang))-(fragLength-posInFrag))>=tolerance)) {
							int newPosition = c - fragStartInSeq - posInFrag - motifPos + queryMotifPos;
							double newScore = -queryScore - 0.0*motifScore;
							if ((newPosition>-1) && (newPosition<energies.length)) {
								if (energies[newPosition]>newScore) {
									energies[newPosition]=newScore;
								}							
							}
						}
					}
				}
			}
		}
	}
	
	public void findingBestEnergy(double[] energies, int seqLen) {
		int[] sortedIndices = AbstractLoopBuilder.findTopMinArray(energies, 3, Double.MAX_VALUE);
		System.out.println("Best energy:" + energies[sortedIndices[0]]);
		System.out.println("In: " + corpus.proteinNames[corpus.protInd[corpus.ungapped[sortedIndices[0]]]] + "  " + 
				corpus.resNum[corpus.ungapped[sortedIndices[0]]]);
		for (int c=0 ; c<seqLen ; c++) {
			System.out.print(corpus.resType[corpus.ungapped[sortedIndices[0]]+c] + " ");
		}
		System.out.println();
	}
	
	private final double[][] confCurvIsite51 = {{10034 , -0.709541 , 5.390276},
			{10040 , -0.644955 , 6.372222},
			{10093 , -1.100516 , 7.867307},
			{10401 , -0.596609 , 1.301971},
			{10523 , -0.777686 , 1.163602},
			{10610 , -0.585547 , 5.298811},
			{10912 , -0.713479 , 4.405523},
			{11007 , -0.708906 , 6.266776},
			{11024 , -0.611499 , 1.347965},
			{11028 , -0.616242 , 1.326293},
			{11040 , -1.125451 , 7.467712},
			{11051 , -1.102956 , 6.347317},
			{11068 , -0.768978 , 5.581519},
			{11070 , -0.349447 , 8.718858},
			{11176 , -0.683162 , 1.386790},
			{11317 , -0.654622 , 1.540378},
			{11386 , -0.249380 , 4.718881},
			{11390 , -0.655639 , 4.781528},
			{11506 , -0.685572 , 4.824707},
			{12001 , -0.503342 , 4.618458},
			{12010 , -0.807776 , 7.461238},
			{12016 , -0.675223 , 6.258407},
			{12018 , -0.658877 , 6.713360},
			{12029 , -0.628176 , 1.523024},
			{12031 , -0.749569 , 6.525108},
			{12046 , -0.521994 , 1.620089},
			{12051 , -0.700398 , 6.521859},
			{12910 , -0.625377 , 4.635984},
			{13001 , -0.470897 , 1.849259},
			{13005 , -0.576061 , 4.832142},
			{13007 , -0.415103 , 5.426950},
			{13008 , -0.643843 , 3.799397},
			{13022 , -0.624481 , 3.678171},
			{13027 , -0.586405 , 1.803015},
			{13066 , -0.686249 , 5.878333},
			{13160 , -0.627469 , 1.761321},
			{13202 , -0.680038 , 5.793248},
			{13907 , -0.669987 , 4.698046},
			{15001 , -0.415451 , 2.273853},
			{15002 , -0.675260 , 4.902541},
			{15004 , -0.424235 , 2.261501},
			{15006 , -0.565805 , 5.380157},
			{15008 , -0.550529 , 3.773436},
			{15015 , -0.498928 , 6.718777},
			{15025 , -0.629965 , 5.124699},
			{15106 , -0.468476 , 2.319948},
			{15286 , -0.555366 , 4.231320},
			{15404 , -0.329752 , 5.230561},
			{15448 , -0.644025 , 4.849470},
			{15548 , -0.392849 , 2.219797},
			{15577 , -0.554382 , 4.210301},
			{15585 , -0.557272 , 5.257100},
			{15610 , -0.626343 , 4.151880},
			{15901 , -0.578161 , 5.134950},
			{15905 , -0.668405 , 4.809428},
			{3009 , -0.666228 , 10.000000},
			{3020 , -0.987108 , 2.430617},
			{3023 , -0.924484 , -0.255711},
			{3025 , -1.188649 , 0.368109},
			{3026 , -1.085394 , 0.361114},
			{3029 , -1.249737 , 0.660052},
			{3030 , -0.958514 , -0.223240},
			{3031 , -1.180219 , -0.538104},
			{3035 , -1.191746 , 0.363862},
			{3041 , -1.178264 , 0.484208},
			{3043 , -1.017365 , -0.368520},
			{3046 , -0.816593 , 4.670798},
			{3047 , -0.533169 , 4.921192},
			{3051 , -0.977285 , -0.274808},
			{3059 , -1.003222 , 0.517793},
			{3060 , -0.692396 , 4.455436},
			{3069 , -0.772572 , 2.718554},
			{3084 , -1.228176 , 0.220872},
			{3108 , -0.932316 , -0.356876},
			{3144 , -0.811247 , 3.007079},
			{3148 , -0.733977 , 4.574355},
			{3162 , -0.854710 , 2.772522},
			{3260 , -0.935011 , 0.843773},
			{4005 , -1.013683 , 0.795590},
			{4008 , -0.884163 , -0.012340},
			{4021 , -0.999161 , 1.221167},
			{4026 , -0.833334 , 0.022639},
			{4031 , -1.151772 , -0.187721},
			{4033 , -0.744366 , 2.552226},
			{4034 , -0.808154 , 3.183349},
			{4058 , -0.958060 , -0.102029},
			{4064 , -0.757390 , 0.119618},
			{4103 , -0.701087 , 3.438346},
			{4114 , -0.908253 , 3.104470},
			{4136 , -0.951736 , 1.074743},
			{4137 , -0.679939 , 3.342503},
			{4143 , -0.872100 , -0.182937},
			{4190 , -0.475749 , 3.380590},
			{4227 , -0.811141 , 1.258533},
			{5032 , -0.707580 , 0.294083},
			{5036 , -0.705085 , 1.988331},
			{5040 , -0.728502 , 0.325518},
			{5050 , -0.730584 , 3.707448},
			{5059 , -1.114219 , 3.595081},
			{5081 , -0.712264 , 3.485656},
			{5097 , -0.869635 , 1.383901},
			{5107 , -0.759990 , 1.876948},
			{5108 , -0.523401 , 3.506878},
			{5116 , -0.308133 , 5.898854},
			{5123 , -0.633450 , 3.499637},
			{5252 , -1.075614 , 0.066026},
			{6037 , -0.762153 , 2.093681},
			{6077 , -0.614443 , 3.258576},
			{6111 , -0.501690 , 4.880779},
			{6139 , -0.688549 , 2.735467},
			{6288 , -0.819138 , 2.554000},
			{6355 , -0.681055 , 2.347458},
			{6513 , -0.830555 , 0.331396},
			{6923 , -0.784016 , 4.619471},
			{7015 , -0.703027 , 3.628658},
			{7037 , -0.762573 , 0.417968},
			{7040 , -0.775994 , 5.119197},
			{7048 , -0.800611 , 5.540289},
			{7135 , -0.686139 , 0.600882},
			{7136 , -0.566722 , 5.718362},
			{7256 , -1.472691 , 4.853507},
			{7312 , -0.600150 , 6.638583},
			{7410 , -0.847794 , 3.976314},
			{7428 , -0.708741 , 3.179843},
			{8103 , -1.377164 , 7.635515},
			{8128 , -0.688313 , 0.783943},
			{8131 , -0.551680 , 3.094829},
			{8165 , -0.479044 , 4.532929},
			{8241 , -0.850323 , 3.557192},
			{8300 , -0.693537 , 4.311580},
			{8332 , -0.753109 , 0.649416},
			{9013 , -0.681057 , 8.684465},
			{9022 , -0.601824 , 1.045877},
			{9024 , -0.755968 , 6.078125},
			{9029 , -0.559283 , 4.934212},
			{9055 , -0.714080 , 5.148932},
			{9073 , -0.623582 , 4.146765},
			{9102 , -0.681481 , 5.506739},
			{9128 , -0.791730 , 6.210047},
			{9164 , -0.729210 , 5.370209},
			{9174 , -0.446279 , 4.246960},
			{9175 , -0.765651 , 3.885276},
			{9480 , -0.627917 , 0.960886},
			{9920 , -0.760585 , 5.492984},
			{9931 , -0.763503 , 7.581025}};

	private final double[][] zScoresIsite51 = {{10034 , -6.847607 , 3.517556},
			{10040 , -6.668476 , 3.232906},
			{10093 , -2.864599 , 1.961345},
			{10401 , -2.609221 , 1.908573},
			{10523 , -2.350284 , 1.699796},
			{10610 , -6.693715 , 3.383707},
			{10912 , -5.218660 , 2.959467},
			{11007 , -6.889216 , 3.298853},
			{11024 , -2.584920 , 1.919996},
			{11028 , -2.600034 , 1.894889},
			{11040 , -2.737991 , 1.918830},
			{11051 , -2.176658 , 1.534762},
			{11068 , -6.296748 , 3.355811},
			{11070 , -6.682289 , 3.327818},
			{11176 , -2.444960 , 1.838117},
			{11317 , -2.528544 , 1.907547},
			{11386 , -2.282340 , 1.899491},
			{11390 , -4.084674 , 2.668808},
			{11506 , -5.952895 , 3.232841},
			{12001 , -6.473442 , 3.439089},
			{12010 , -6.572576 , 3.459305},
			{12016 , -5.871665 , 3.228496},
			{12018 , -6.540392 , 3.235332},
			{12029 , -2.627359 , 1.954237},
			{12031 , -6.592436 , 3.371792},
			{12046 , -2.852956 , 2.106656},
			{12051 , -6.882414 , 3.366908},
			{12910 , -5.733478 , 3.197044},
			{13001 , -3.038994 , 2.258178},
			{13005 , -5.360933 , 3.128956},
			{13007 , -7.165519 , 3.763979},
			{13008 , -4.590991 , 2.899541},
			{13022 , -6.140065 , 3.417717},
			{13027 , -2.778989 , 2.081679},
			{13066 , -4.236603 , 2.748714},
			{13160 , -2.734225 , 2.055591},
			{13202 , -7.130162 , 3.499750},
			{13907 , -5.811207 , 3.240961},
			{15001 , -3.782053 , 2.669762},
			{15002 , -6.483142 , 3.465139},
			{15004 , -3.653853 , 2.604583},
			{15006 , -8.276078 , 3.977731},
			{15008 , -5.509302 , 3.257752},
			{15015 , -6.013135 , 3.471785},
			{15025 , -6.526935 , 3.484634},
			{15106 , -3.512384 , 2.542615},
			{15286 , -6.085565 , 3.403910},
			{15404 , -3.053433 , 2.408547},
			{15448 , -6.570263 , 3.520306},
			{15548 , -3.948667 , 2.751239},
			{15577 , -6.547284 , 3.547934},
			{15585 , -8.363054 , 4.014182},
			{15610 , -6.299733 , 3.447692},
			{15901 , -8.615384 , 4.077225},
			{15905 , -6.130119 , 3.375506},
			{3009 , -0.996313 , 0.216743},
			{3020 , -2.104262 , 1.413926},
			{3023 , -1.644246 , 0.990760},
			{3025 , -1.231337 , 0.563812},
			{3026 , -1.390032 , 0.728801},
			{3029 , -1.235353 , 0.596664},
			{3030 , -1.608773 , 0.959048},
			{3031 , -1.554811 , 0.891147},
			{3035 , -1.246342 , 0.579858},
			{3041 , -1.341009 , 0.698574},
			{3043 , -1.621182 , 0.950100},
			{3046 , -2.728948 , 1.676560},
			{3047 , -3.079673 , 1.735260},
			{3051 , -1.551402 , 0.884891},
			{3059 , -1.326359 , 0.690075},
			{3060 , -1.824226 , 1.280327},
			{3069 , -4.002210 , 2.280605},
			{3084 , -1.243183 , 0.563298},
			{3108 , -1.663460 , 0.977153},
			{3144 , -3.377460 , 1.999909},
			{3148 , -3.197028 , 1.876681},
			{3162 , -1.293903 , 0.793994},
			{3260 , -1.324962 , 0.722463},
			{4005 , -1.503943 , 0.877671},
			{4008 , -1.760707 , 1.108760},
			{4021 , -1.392339 , 0.824138},
			{4026 , -1.954635 , 1.267215},
			{4031 , -1.787380 , 1.133422},
			{4033 , -3.217368 , 2.021857},
			{4034 , -4.107077 , 2.359662},
			{4058 , -1.755419 , 1.099582},
			{4064 , -1.989098 , 1.299391},
			{4103 , -3.245176 , 1.904875},
			{4114 , -3.704694 , 2.215693},
			{4136 , -1.497383 , 0.928071},
			{4137 , -2.833338 , 1.895349},
			{4143 , -1.903033 , 1.194301},
			{4190 , -1.522811 , 1.054536},
			{4227 , -1.611718 , 1.030244},
			{5032 , -2.107039 , 1.414500},
			{5036 , -1.713191 , 1.170974},
			{5040 , -1.953403 , 1.309309},
			{5050 , -3.863513 , 2.257510},
			{5059 , -1.808235 , 1.286009},
			{5081 , -3.951083 , 2.343660},
			{5097 , -1.803017 , 1.193491},
			{5107 , -1.663037 , 1.134472},
			{5108 , -1.574772 , 1.131101},
			{5116 , -1.387678 , 0.951964},
			{5123 , -3.626397 , 2.287832},
			{5252 , -1.798925 , 1.165699},
			{6037 , -1.807129 , 1.272091},
			{6077 , -2.038555 , 1.460936},
			{6111 , -4.330777 , 2.459275},
			{6139 , -1.789474 , 1.299613},
			{6288 , -1.677890 , 1.181371},
			{6355 , -2.001512 , 1.423557},
			{6513 , -2.178990 , 1.478817},
			{6923 , -5.108826 , 2.664107},
			{7015 , -4.246562 , 2.589651},
			{7037 , -2.288928 , 1.557735},
			{7040 , -4.895546 , 2.695852},
			{7048 , -3.337810 , 2.202647},
			{7135 , -2.409822 , 1.674053},
			{7136 , -4.896943 , 2.873002},
			{7256 , -1.559533 , 1.042531},
			{7312 , -2.882852 , 2.063392},
			{7410 , -3.364272 , 2.163098},
			{7428 , -1.695173 , 1.241572},
			{8103 , -1.631408 , 1.190696},
			{8128 , -2.444166 , 1.717216},
			{8131 , -2.803965 , 1.926231},
			{8165 , -6.054040 , 3.152139},
			{8241 , -3.776814 , 2.403194},
			{8300 , -3.996215 , 2.455778},
			{8332 , -2.156274 , 1.510102},
			{9013 , -5.362850 , 3.047867},
			{9022 , -2.434746 , 1.761392},
			{9024 , -6.008761 , 3.014731},
			{9029 , -5.500932 , 3.086649},
			{9055 , -5.760630 , 2.928271},
			{9073 , -6.046406 , 3.249904},
			{9102 , -6.395461 , 3.187247},
			{9128 , -5.794061 , 2.933072},
			{9164 , -5.729203 , 3.017672},
			{9174 , -3.432478 , 2.274681},
			{9175 , -4.167425 , 2.605789},
			{9480 , -2.542178 , 1.807700},
			{9920 , -5.712082 , 2.997844},
			{9931 , -6.423223 , 3.159671}};
	
	public double[] phiPsiInSeq(int res) {
		return phiPsiInSeq[res];
	}
	
	public static void main(String[] args) {
        
//		Command command = commands.firstWord(PARAMETERS_DIRECTORY);
//        String parametersPath = command.secondWord();
//        command = commands.firstWord(ROTAMER_LIBRARY);
//        parametersFileName = parametersFileName+"/"+command.secondWord();

		IsiteLib lib = new IsiteLib(args[0]);
		System.out.println("Loaded Isite");
		Corpus corpus = new Corpus(args[1]);
		System.out.println("Loaded corpus");
		lib.setCorpus(corpus);
//		lib.calcMotifsInCorpus();
		lib.loadMotifInCorpus("motifAssignment.dat");
		Protein prot = new Protein("C:\\Users\\Nir\\Loop_Building_Project\\Pisces_REDUCE\\1EDQA.pdb",
				new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		lib.loadAndAnalyzeSequenceIntoLib(Protein.getSeqOfProt(prot, 0, prot.lastResidue()));
		for (int pos=0 ; pos<prot.lastResidue() ; pos++) {
			lib.bestMotifInSeqPos(pos);
		}
//		int[] querySeq = {17,4,8,12,15,16,5,2,9,7};
//		corpus.buildUngappedArray(querySeq.length);
//		double[] energies = new double[corpus.ungapped.length];
//		lib.setThreadingEnergies(querySeq, 1, 6, energies);
//		lib.findingBestEnergy(energies, querySeq.length);
//		for (int c=0 ; c<querySeq.length ; c++) {
//			System.out.print(querySeq[c] + " ");
//		}
//		System.out.println();
		

		
//		lib.calcMotifsInCorpus();
//		lib.findCenterOfMotif(79);
//		for (int c=0 ; c<lib.size() ; c++) {
//			lib.calcSuccess(c);
//		}
	}
	
	
	
}
