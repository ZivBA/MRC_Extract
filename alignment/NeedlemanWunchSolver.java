package alignment;

public class NeedlemanWunchSolver {
	
	protected Sequence seq1;
	protected Sequence seq2;
	protected ScoringScheme scoring;
	private double[][] Smatrix;
	private int[][] Dmatrix;
	private double qualityIndex;
	protected AlignmentReturnStructure bestAlign;
	
	
	public  NeedlemanWunchSolver(Sequence seq1, Sequence seq2, ScoringScheme scoring) {
		this.scoring = scoring;
		this.seq1 = seq1;
		this.seq2 = seq2;
		compute();
	}
	
	public void printSmatrix() {
		for (int i=0; i<Smatrix.length ; i++) {
			for (int j=0; j<Smatrix[i].length ; j++) {
				System.out.print(Math.round(Smatrix[i][j]*100)/100.0 + " ");
			}
			System.out.println();
		}
	}
	
	public void printDmatrix() {
		for (int i=0; i<Dmatrix.length ; i++) {
			for (int j=0; j<Dmatrix[i].length ; j++) {
				System.out.print(Math.round(Dmatrix[i][j]*100)/100.0 + " ");
			}
			System.out.println();
		}
	}
	
	public void printAlignment() {
		System.out.println(logAlignment());
	}
	String logAlignment() {
		bestAlign = backtrack(Smatrix.length-1, Smatrix[0].length-1);
		return("Score: " + alignmentScore() + "  Matches: " + bestAlign.matches + "\n" + bestAlign.S1 + "\n" + bestAlign.S2 + "\n" + "Quality " +
				"Index: " + qualityIndex);
	}
	
	public double alignmentScore() {
		return Smatrix[Smatrix.length-1][Smatrix[0].length-1];
	}
	
	private void compute() {
		Smatrix = new double[seq2.length()+1][seq1.length()+1];
		Dmatrix = new int[seq2.length()+1][seq1.length()+1];
		double[] RowMaxPenaltyScore = new double[seq1.length()+1];
		double[] ColMaxPenaltyScore = new double[seq2.length()+1];
		int[] RowWhereBestPenaltyScore = new int[seq1.length()+1];
		int[] ColWhereBestPenaltyScore = new int[seq2.length()+1];
		for (int c=0; c<RowMaxPenaltyScore.length ; c++) {
			Smatrix[0][c] = -10000.0*c; // A large number because initially no gaps are allowed on the Rval side.
			RowMaxPenaltyScore[c] = Double.NEGATIVE_INFINITY;
			Dmatrix[0][c] = -1;
		}
		for (int c=0; c<ColMaxPenaltyScore.length ; c++) {
			Smatrix[c][0] = -0.25*c;
			ColMaxPenaltyScore[c] = Double.NEGATIVE_INFINITY;
			Dmatrix[c][0] = 1;
		}
		
		// Computing the rest of the matrix
		for (int i=1; i<=seq2.length() ; i++) {
			for (int j=1 ; j<=seq1.length() ; j++) {
				double diagScore = Math.max(scoring.score(seq1.get(j-1) , seq2.get(i-1)) ,
						Math.max(seq1.get(j-1).gapOpeningScore()+seq1.get(j-1).gapAligningScore() , seq2.get(i-1).gapOpeningScore()+seq2.get(i-1).gapAligningScore())); // Positive gap-penalty correction
				Smatrix[i][j]=Smatrix[i-1][j-1]+diagScore;
				Dmatrix[i][j]=0;
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
				// This part is only possible in positive gap-penalty option.
				if (Dmatrix[i][j]==0) {
					Smatrix[i][j]=Smatrix[i-1][j-1]+scoring.score(seq1.get(j-1) , seq2.get(i-1));
				}
				// End of positive gap-penalty correction
			} // of j
		}  // of i
		
		qualityIndex = Smatrix[Smatrix.length-1][Smatrix[0].length-1] / seq1.size();
	} // Of compute
	
	
	/**
	 * Backtracking from a certain position in the matrix.
	 **/
	private AlignmentReturnStructure backtrack(int startI, int startJ) {
		String S1 = "";
		String S2 = "";
		while ((startI>0) || (startJ>0)) {
			if (Dmatrix[startI][startJ]==0) {
				S1=seq1.get(startJ-1).string() + S1;
				S2=seq2.get(startI-1).string() + S2;
				startI--;
				startJ--;
			}
			else if (Dmatrix[startI][startJ]<0) {
				for (int c=startJ ; c>startJ-Math.abs(Dmatrix[startI][startJ]) ; c--) {
					S1=seq1.get(c-1).string() + S1;
					S2=seq2.get(Math.max(0 , startI-1)).gapString() + S2;  // gaps in beginning are as defined for first position
				}
				startJ-=Math.abs(Dmatrix[startI][startJ]);
			}
			else {
				for (int c=startI ; c>startI-Math.abs(Dmatrix[startI][startJ]) ; c--) {
					S1=seq1.get(Math.max(0 , startJ-1)).gapString() + S1;  // gaps in beginning are as defined for first position
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
	
	public double getQualityIndex() {
		return qualityIndex;
	}
	
	public class AlignmentReturnStructure {
		public int matches;
		public String S1;
		public String S2;
	}
	
}
