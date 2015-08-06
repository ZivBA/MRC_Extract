package meshi.applications.loopBuilding;


public class BasicLoopResult {
		protected double evEnergy;
		protected double propEnergy;
		protected double bbHBEnergy;
		protected double score;		
		protected double rms;		
		protected double[][] coors = null; // This pointer refer to the output of the "saveLoopCoordinates()" method
		protected double[][] myPhiPsi = null; // This pointer refer to the output of the "savePhiPsiOfLoop()" method
		protected int[] fragRank = null; // From where in the library did the fragments which makes the solution come?
		private double[][] bbCoors = null; // for N residues this is a [3N][3] array. Indexing is N,Ca,C,N,Ca,C... and the other index is {x,y,z}

		public BasicLoopResult() {};
		
		public BasicLoopResult(double evEnergy, double propEnergy, double bbHBEnergy, double score, double rms, 
				int[] fragRank_p, double[][] coors, double[][] bbCoors, double[][] myPhiPsi) {
			this.evEnergy = evEnergy;
			this.propEnergy = propEnergy;
			this.bbHBEnergy = bbHBEnergy;
			this.score = score;
			this.rms = rms;
			this.coors = coors;
			this.bbCoors = bbCoors;
			this.myPhiPsi = myPhiPsi;
			fragRank = new int[fragRank_p.length];
			for (int c=0 ; c<fragRank_p.length ; c++)
				fragRank[c] = fragRank_p[c];
		}

		public double[][] getBBcoors() {return bbCoors;}
		
		public double calcRMS(BasicLoopResult other) {
			return calcRMS(other.bbCoors);
		}

		/**
		 * The real RMS value
		 */
		public double calcRMS(double[][] other) {
			if (bbCoors[0].length != other[0].length)
				throw new RuntimeException("ERROR: the bbcoors arrays in both instances are not the same length.");
			double totRms = 0.0;
			for (int c=0 ; c<bbCoors[0].length ; c++) 
				totRms += 
					((bbCoors[0][c]-other[0][c])*(bbCoors[0][c]-other[0][c]) +
					(bbCoors[1][c]-other[1][c])*(bbCoors[1][c]-other[1][c]) +
					(bbCoors[2][c]-other[2][c])*(bbCoors[2][c]-other[2][c]));
			return Math.sqrt(totRms/bbCoors[0].length);
		}
		
		/**
		 * The maximal deviation
		 */
//		public double calcRMS(double[][] other) {
//			if (bbCoors[0].length != other[0].length)
//				throw new RuntimeException("ERROR: the bbcoors arrays in both instances are not the same length.");
//			double maxDev = 0.0;
//			double dev = 0.0;
//			for (int c=0 ; c<bbCoors[0].length ; c++) { 
//				dev = ((bbCoors[0][c]-other[0][c])*(bbCoors[0][c]-other[0][c]) +
//					(bbCoors[1][c]-other[1][c])*(bbCoors[1][c]-other[1][c]) +
//					(bbCoors[2][c]-other[2][c])*(bbCoors[2][c]-other[2][c]));
//				if (dev>maxDev)
//					maxDev = dev;
//			}
//			return Math.sqrt(maxDev);
//		}
		
		public void setScore(double newScore) {
			score = newScore;
		}
		
}
