package alignment;

import java.util.Arrays;

public class SmithWatermanSolver {

	protected Sequence seq1;
	protected Sequence seq2;
	protected ScoringScheme scoring;
	private double[][] Smatrix;
	private int[][] Dmatrix;
	
	public  SmithWatermanSolver(Sequence seq1, Sequence seq2, ScoringScheme scoring) {
		this.scoring = scoring;
		this.seq1 = seq1;
		this.seq2 = seq2;
		compute();
	}

	
	/**
	 * Returns the top N scores 
	 **/
	public double[] topScores(int N) {
		double[] scores = new double[seq1.size()*seq2.size()];
		int c=0;
		for (int i=1 ; i<Smatrix.length ; i++) {
			for (int j=1 ; j<Smatrix[0].length ; j++) {
				scores[c] = Smatrix[i][j];
				c++;
			}
		}
		Arrays.sort(scores);
		double[] out = new double[Math.min(N, scores.length)];
		for (c=0 ; c<out.length ; c++) {
			out[c] = scores[scores.length-c-1];
		}
		return out;
	}

	
	/**
	 * Finds the alignment with the highest number of matches in the top N scores.  
	 **/
	public void mostMatchesAlignment(int N) {
		int bestMatch = -1;
		int bestMatchI = -1;
		int bestMatchJ = -1;
		double[] topScores = topScores(N);
		for (int i=1 ; i<Smatrix.length ; i++) {
			for (int j=1 ; j<Smatrix[0].length ; j++) {
				if (Smatrix[i][j]>=topScores[topScores.length-1]) {
					AlignmentReturnStructure alignment = backtrack(i, j);
					if (alignment.matches>bestMatch) {
						bestMatch = alignment.matches;
						bestMatchI = i;
						bestMatchJ = j;
					}
				}
			}
		}
		for (int c=0 ; c<topScores.length ; c++) {
			if (Smatrix[bestMatchI][bestMatchJ]>topScores[c]) {
				AlignmentReturnStructure bestMatchAlignment = backtrack(bestMatchI, bestMatchJ);
				System.out.println("Best match - rank: " + c + "  Score: " + topScores[c-1] + "  (" + topScores[topScores.length-1] + ")  " + "  Matches: " + bestMatchAlignment.matches + 
						"\n" + bestMatchAlignment.S1 + "\n" + bestMatchAlignment.S2 + "\n");
				return;
			}
		}
	}

	public void printBestAlignment() {
		double bestScore = -1;
		int bestI = -1;
		int bestJ = -1;
		for (int i=0 ; i<Smatrix.length ; i++) {
			for (int j=0 ; j<Smatrix[0].length ; j++) {
				if (Smatrix[i][j]>bestScore) {
					bestScore = Smatrix[i][j];
					bestI = i;
					bestJ = j;
				}
			}
		}
		AlignmentReturnStructure bestAlign = backtrack(bestI, bestJ);
		System.out.println("Best score: " + bestScore + "  Matches: " + bestAlign.matches + "\n" + bestAlign.S1 + "\n" + bestAlign.S2);
	}

	public double bestAlignmentScore() {
		double bestScore = -1;
		for (int i=0 ; i<Smatrix.length ; i++) {
			for (int j=0 ; j<Smatrix[0].length ; j++) {
				if (Smatrix[i][j]>bestScore) {
					bestScore = Smatrix[i][j];
				}
			}
		}
		return bestScore;
	}

	private void compute() {
		Smatrix = new double[seq2.length()+1][seq1.length()+1];
		Dmatrix = new int[seq2.length()+1][seq1.length()+1];
		double[] RowMaxPenaltyScore = new double[seq1.length()+1];
		double[] ColMaxPenaltyScore = new double[seq2.length()+1];
		int[] RowWhereBestPenaltyScore = new int[seq1.length()+1];
		int[] ColWhereBestPenaltyScore = new int[seq2.length()+1];
		for (int c=1; c<RowMaxPenaltyScore.length ; c++) {
			RowMaxPenaltyScore[c] = Double.NEGATIVE_INFINITY;
			Dmatrix[0][c] = -1;
		}
		for (int c=1; c<ColMaxPenaltyScore.length ; c++) {
			ColMaxPenaltyScore[c] = Double.NEGATIVE_INFINITY;
			Dmatrix[c][0] = 1;
		}

		// Computing the rest of the matrix
		for (int i=1; i<=seq2.length() ; i++) {
			for (int j=1 ; j<=seq1.length() ; j++) {
				// Match
				double diagScore = scoring.score(seq1.get(j-1) , seq2.get(i-1));
				if (Smatrix[i-1][j-1]+diagScore>=Smatrix[i][j]) {
					Smatrix[i][j]=Smatrix[i-1][j-1]+diagScore;
					Dmatrix[i][j]=0;					
				}
				// gap along the row
				double extendedGapScore = ColMaxPenaltyScore[i]+seq2.get(i-1).gapAligningScore();
				double newGapScore = Smatrix[i][j-1]+seq2.get(i-1).gapOpeningScore()+seq2.get(i-1).gapAligningScore();
				if (Dmatrix[i][j-1] < 0) { // This is neccessary when there is a positive gap-penalty (an incentive for gap opening)   - No NEW gap OPENING in an open gap
					newGapScore = Double.NEGATIVE_INFINITY;
				}
				if (newGapScore>extendedGapScore) { // Better to open a new gap here then extended the best so far
					if (newGapScore>Smatrix[i][j]) {
						Smatrix[i][j]=newGapScore;
						Dmatrix[i][j]=-1; // Going along row is marked with negative numbers
						ColMaxPenaltyScore[i] = newGapScore;
						ColWhereBestPenaltyScore[i] = j;
					}
				}
				else { // Better to extend an existing gap
					if (extendedGapScore>Smatrix[i][j]) {
						Smatrix[i][j]=extendedGapScore;
						Dmatrix[i][j]=-(j-ColWhereBestPenaltyScore[i] +1); // Going along row is marked with negative numbers
						ColMaxPenaltyScore[i] = extendedGapScore;
					}
				}
				// gap along the Column
				extendedGapScore = RowMaxPenaltyScore[j]+seq1.get(j-1).gapAligningScore();
				newGapScore = Smatrix[i-1][j]+seq1.get(j-1).gapOpeningScore()+seq1.get(j-1).gapAligningScore();
				if (Dmatrix[i-1][j] > 0) { // This is neccessary when there is a positive gap-penalty (an incentive for gap opening)   - No NEW gap OPENING in an open gap
					newGapScore = Double.NEGATIVE_INFINITY;
				}
				if (newGapScore>extendedGapScore) { // Better to open a new gap here then extended the best so far
					if (newGapScore>Smatrix[i][j]) {
						Smatrix[i][j]=newGapScore;
						Dmatrix[i][j]=1; // Going along column is marked with positive numbers
						RowMaxPenaltyScore[j] = newGapScore;
						RowWhereBestPenaltyScore[j] = i;
					}
				}
				else { // Better to extend an existing gap
					if (extendedGapScore>Smatrix[i][j]) {
						Smatrix[i][j]=extendedGapScore;
						Dmatrix[i][j]=(i-RowWhereBestPenaltyScore[j] +1); // Going along column is marked with positive numbers
						RowMaxPenaltyScore[j] = extendedGapScore;
					}
				}
			} // of j
		}  // of i		
	} // Of compute

	
	/**
	 * Backtracking from a certain position in the matrix.
	 **/
	private AlignmentReturnStructure backtrack(int startI, int startJ) {
		String S1 = "";
		String S2 = "";
		while ((startI>0) && (startJ>0) && (Smatrix[startI][startJ]>0)) {
			if (Dmatrix[startI][startJ]==0) {
				S1=seq1.get(startJ-1).string() + S1;
				S2=seq2.get(startI-1).string() + S2;
				startI--;
				startJ--;
			}
			else if (Dmatrix[startI][startJ]<0) {
				for (int c=startJ ; c>startJ-Math.abs(Dmatrix[startI][startJ]) ; c--) {
					S1=seq1.get(c-1).string() + S1;
					S2=seq2.get(startI-1).gapString() + S2;
				}
				startJ-=Math.abs(Dmatrix[startI][startJ]);
			}
			else {
				for (int c=startI ; c>startI-Math.abs(Dmatrix[startI][startJ]) ; c--) {
					S1=seq1.get(startJ-1).gapString() + S1;
					S2=seq2.get(c-1).string() + S2;
				}
				startI-=Math.abs(Dmatrix[startI][startJ]);
			}
		}
		
		// Formating the return structure
		AlignmentReturnStructure out = new AlignmentReturnStructure();
		int matches=0;
		for (int c=0; c<S1.length() ; c++) {
			if (S1.charAt(c)==S2.charAt(c)) {
				matches++;
			}
		}
		out.S1 = S1;
		out.S2 = S2;
		out.matches = matches;
		return out;
	}

	public class AlignmentReturnStructure {
		public int matches;
		public String S1;
		public String S2;
	}
	
}

