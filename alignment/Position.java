package alignment;

public interface Position {
	
	/*
	 * The penalty involving a gap opening against this position
	 */
	public double gapOpeningScore();
	
	/*
	 * The penalty involving aligning this position to a gap
	 */
	public double gapAligningScore();
	
	/*
	 * A string for printing out the position in the alignment
	 */
	public String string();

	/*
	 * A string for printing out a gap
	 */
	public String gapString();

	
}
