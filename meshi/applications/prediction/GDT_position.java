package meshi.applications.prediction;

import meshi.util.overlap.Overlap;

public class GDT_position {
	private double[][] C1;     
	private double[][] C2;  
	private double[][] copyOfInitialC2,copyOfInitialC1;
	public int[] subset;
	private double cutoff;

	public GDT_position(int sub, int place, double max_distance, double[][] C1, double[][] C2)  
	{
		this.C1 = C1;
		this.C2 = C2;
		// Making hard copies of C2 and C1, for cases where the RMS calculations are not successful and garbage the arrays.
		copyOfInitialC2 = new double[3][C2[0].length];
		copyOfInitialC1 = new double[3][C1[0].length];
		for (int c=0 ; c<C2[0].length ; c++) {
			copyOfInitialC2[0][c] = C2[0][c];
			copyOfInitialC2[1][c] = C2[1][c];
			copyOfInitialC2[2][c] = C2[2][c];
			copyOfInitialC1[0][c] = C1[0][c];
			copyOfInitialC1[1][c] = C1[1][c];
			copyOfInitialC1[2][c] = C1[2][c];
		}
		subset = new int[sub];
		for (int i=0 ; i<sub ; i++) {
			subset[i] = place + i;
		}
		cutoff = max_distance*max_distance;
	}

	public void find_best_conformation(int num_of_loops)
	{
		int[] oldset;
		boolean conv;
		int loopCount=0;
		do {
			conv = true;
			oldset = subset;
			try {
				doIter();
			}
			catch (Exception e) {
				System.out.println("Warning: a problem was encountered during the GDT search." + 
						" If this problem occurred only once or twice for this model, then the final GDT result" +
						" is reliable. However if this warning occured frequently on this model, then there is "+
						"some ill-defined properties to it (like all atoms are co-planar, for example). In this case " +
						"the GDT result is not reliable. The cause of warning is: \n\n\n" + e.getMessage() +
						"\n\n\n");
				// Because some RMS errors tends to garbage C2, we return it its original values.
				for (int c=0 ; c<C2[0].length ; c++) {
					C2[0][c] = copyOfInitialC2[0][c];
					C2[1][c] = copyOfInitialC2[1][c];
					C2[2][c] = copyOfInitialC2[2][c];
					C1[0][c] = copyOfInitialC1[0][c];
					C1[1][c] = copyOfInitialC1[1][c];
					C1[2][c] = copyOfInitialC1[2][c];
				}
			}
			if (oldset.length != subset.length)
				conv = false;
			else {
				for (int j=0 ; (j<subset.length)&&conv ; j++) {
					if (subset[j] != oldset[j])
						conv = false;
				}
			}      
			if (subset.length < 3)
				conv = true; 
			if (oldset.length > subset.length) {
				subset = oldset;
				conv = true;
			}
			loopCount++;
		} while ((loopCount<num_of_loops) && (!conv));
	}

	//----------------------------------------------------
	private void doIter()
	{
		double rms;
		rms = Overlap.rmsPartial(C1,C2,subset);
		if ((rms <0) || Double.isNaN(rms))
			throw new RuntimeException("\n\nincorrect rms: " + rms + "\n\n");
		updateSubset();
	}

	//----------------------------------------------------
	/**Updateing the subset**/
	private void updateSubset()
	{
		int[] tmpset = new int[C1[0].length];  
		int i,j;
		j = 0;
		for (i=0 ; i<C1[0].length ; i++){
			if (((C1[0][i] - C2[0][i])*(C1[0][i] - C2[0][i])+
					(C1[1][i] - C2[1][i])*(C1[1][i] - C2[1][i])+
					(C1[2][i] - C2[2][i])*(C1[2][i] - C2[2][i])) < cutoff) {
				tmpset[j] = i;
				j++;
			}
		}
		subset = new int[j];
		for (i=0 ; i<j ; i++){
			subset[i] = tmpset[i];
		}
	}

} //class GDT_position