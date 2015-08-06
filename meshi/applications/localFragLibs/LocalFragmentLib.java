package meshi.applications.localFragLibs;

import java.text.DecimalFormat;
import java.util.Arrays;

import meshi.applications.loopBuilding.AbstractLoopBuilder;
import meshi.applications.corpus.Corpus;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.geometry.ResidueBuilder;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.parameters.AtomTypes;
import meshi.parameters.MeshiPotential;
import meshi.parameters.Residues;
import meshi.util.KeyWords;
import meshi.util.MeshiProgram;
import meshi.util.Isite.IsiteLib;

/**
 *The purpose of this class is is class is to build a fragment library that is compatible with a certain 
 *sequence. The lib is applicable to a corpus instance. The sequence is given to the constructor as an array
 *of integers (in the range: 0-19). The constructor creates a library with at most libSize fragments that are 
 *different from each other by at least a certain RMS curoff. The final size of the lib could be smaller than
 *the requested if this cutoff is stringent. In general the RMS (for fragment diversity purposes) is calculated 
 *not on the entire sequence, but rather on an inner subsequence. I found out that when you do threading
 *on a long sequence the results are more accurate for the inner subsequences. The ranking of fragments in the lib 
 *is currently by propensity energy and pair-wise sequence alignment using BLOSUM62. 
 **/

public class LocalFragmentLib implements Residues,AtomTypes,CompositeTorsionsDefinitions,MeshiPotential,
KeyWords { 

	// Constants
	private static final double W_PROP_ENERGY =  1.5; 
	private static final double W_BLOSUM62_ENERGY = -1.0; // Note the minus sign!
	
	private Corpus corpus = null;
	private int fragL = 0;
	private int libSize = 0;
	private int desiredLibSize = -1;  // The desired lib size set by the user. 
	private double rmsCutOff;
	private int[] seq = null;
	private int[] libOrig = null;  // This is the lib core. The indices here point to the fragments in the corpus. The size of this array is the lib size.
	private double[] libEnergy = null;  // This is the energy ascociated with a fragment. The size of this array is the lib size.
	protected boolean[] similarToLib = null; // This marks (true values) fragments that are already selected, or are very similar to fragments alreay in the lib. The size of this array is as curpus.ungapped
	private int trueFragStart = -999;  // The place in the truncated (no extension) sequence where the true frag (no overlap) starts
	private int trueFragEnd = -999;  // The place in the truncated (no extension) sequence where the true frag (no overlap) ends
	private int fragStartInSeq = -999;	// See the definition of this field in the constructor documentation.
	private int fragEndInSeq = -999;	// See the definition of this field in the constructor documentation.
	private double[] threadingEnergies = null; // The propensity+BLOSUM62 energies of each fragment in the ungapped corpus.
	private int statusPG = 0;  // 3-no PG ; 2-stringent PG only on the fragment; 1-full stringent PG 
	private LocalFragmentLib prevLib  = null;

//	These fields are needed when the library is constructed so that it is compatible with the stalks.
	private Protein prot = null;
	private int residueFragStart = -999;
	private int overlap = -999;
	private int manner = -999; 
	protected double overlapSimilarityTH = 0.0;
	
	/* This constructor is able to construct a frag-lib, with some lateral addendums on the fragment as follows:
	 * 
	 * |-------------------------------------------SEQ------------------------------------------------------------------------|
	 * [   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ]
	 *                          |--overlap----||--------------------frag--------------||--overlap----|
	 * |---------expandThreadingN-------------|                                       |-----------------expandThreadingC------|
	 * 
	 * Where "frag" is what we are interested in. We make sure, all the "frag" parts are at least 'rmsCutoff' mutually different than 
	 * all the other "frags" in the library. In this constructor this difference is calculated after the OVERLAP regions 
	 * are superimposed. Overlap may be on one side of "frag" or both, depending on 'manner' (which can be -1,0,1). The overlap
	 * must also have an RMS of less than 'overlapSimilarityTH' to the protein in the overlapping region.
	 * 
	 * The parameter 'fragStartInSeq' is the index in 'seq' (begining of seq is 0) where 'frag' or 'overlap' (the first of the two) begin.
	 * The parameter 'fragEndInSeq' is the index in 'seq' (begining of seq is 0) where 'frag' or 'overlap' (the latter of the two) end.
	 * The parameter 'residueFragStart' should be the residue number in 'prot' where 'frag' or 'overlap' (the first of the two) begin.
	 * 
	 */
	public LocalFragmentLib(Corpus corpus , int[] seq, int  desiredLibSize,  double rmsCutOff , int fragStartInSeq , int fragEndInSeq , 
			Protein prot, int residueFragStart, int overlap, int manner, double overlapSimilarityTH, int statusPG) {
		this.corpus = corpus;
		this.seq = seq;
		this.desiredLibSize = desiredLibSize;
		this.rmsCutOff = rmsCutOff;
		this.prot = prot;
		this.residueFragStart = residueFragStart;
		this.overlap = overlap;
		this.manner = manner;
		this.overlapSimilarityTH = overlapSimilarityTH;
		this.fragStartInSeq = fragStartInSeq;
		this.fragEndInSeq = fragEndInSeq;
		this.statusPG = statusPG;

		// Where is the frag in the final sequence? i.e. not considering the overlap(s). This part is required for inheriting classes. It must come before build lib.
		if (manner==-1) {  // N term
			trueFragStart = overlap;
			trueFragEnd = fragL-1;
		}
		else if (manner==1) {  // C term
			trueFragStart = 0;
			trueFragEnd = fragL-1-overlap;
		}
		else if (manner==0) {  // both sides
			trueFragStart = overlap;
			trueFragEnd = fragL-1-overlap;
		}

		//  Calculating the threading energy according to the expanded sequence
		corpus.buildUngappedArray(seq.length);
		threadingEnergies = calculateThreadingEnergy(); 
		
		// After the threading energy was calculated, the sequence can be cut to the essentials. 
		this.seq = new int[fragEndInSeq-fragStartInSeq+1];
		fragL = this.seq.length;
		for (int c=0 ; c<fragL ; c++)
			this.seq[c] = seq[fragStartInSeq+c];
		
		// And now building the lib
		System.out.println("Start building lib.");
		initializeSimilarToLib();
		// I added this on 10/1/2009 to prevent bad PG placements
		for (int c=0 ; c<similarToLib.length ; c++)
			if (Double.isNaN(threadingEnergies[c]))
				similarToLib[c] = true;
		buildLib(true);
		similarToLib = null;  // We don't need it anymore, and it takes quite a lot of memory
		threadingEnergies = null; // We don't need it anymore, and it takes quite a lot of memory		
	}	
		
	

	/**
	 * This constructor is to be used when the conformation on which this frag is to be superimposed is not known during construction. Nontheless 
	 * there should still be defined 'overlap' and 'manner' to this lib, even though they are not explicitly given. 'overlap' and 'manner' come into play in the 
	 * definitions of 'fragStartInSeq' and 'fragEndInSeq' which should be set as for the previous constructor. 
	 **/
	public LocalFragmentLib(Corpus corpus , int[] seq, int desiredLibSize , double rmsCutOff , int fragStartInSeq, int fragEndInSeq, 
			int manner, int overlap, int statusPG, double overlapSimilarityTH, LocalFragmentLib prevLib) {
		this.corpus = corpus;
		this.seq = seq;
		this.fragL = seq.length;
		this.desiredLibSize = desiredLibSize;
		this.rmsCutOff = rmsCutOff;
		this.fragStartInSeq = fragStartInSeq;
		this.fragEndInSeq = fragEndInSeq;
		this.overlap = overlap;
		this.manner = manner;
		this.statusPG = statusPG;
		this.prevLib = prevLib;
		this.overlapSimilarityTH = overlapSimilarityTH;
		
		//  Calculating the threading energy according to the expanded sequence
		corpus.buildUngappedArray(seq.length);
		threadingEnergies = calculateThreadingEnergy();

		// After the threading energy was calculated, the sequence can be cut to the essentials. 
		this.seq = new int[fragEndInSeq-fragStartInSeq+1];
		fragL = this.seq.length;
		for (int c=0 ; c<fragL ; c++)
			this.seq[c] = seq[fragStartInSeq+c];

		// And now building the lib
		System.out.println("Start building lib.");
		initializeSimilarToLib();
		// I added this on 10/1/2009 to prevent bad PG placements
		for (int c=0 ; c<similarToLib.length ; c++)
			if (Double.isNaN(threadingEnergies[c]))
				similarToLib[c] = true;
		buildLib(false);
		similarToLib = null;  // We don't need it anymore, and it takes quite a lot of memory
		threadingEnergies = null; // We don't need it anymore, and it takes quite a lot of memory
	}	
	
	
	// This calculate the threading energy, according to the sequence and the the 'corpus.ungapped' array (built to the sequence length) 
	private double[] calculateThreadingEnergy() {		
		double[] energies = new double[corpus.ungapped.length];
				
		// Checking the PG status
		if ((statusPG<1) | (statusPG>3))
			throw new RuntimeException("The PG status must be {1,2 or 3}.");

		// Finding pre-Prolines in the sequence
		boolean[] pp = new boolean[seq.length];
		for (int e=0; e<seq.length-1 ; e++)
			if ((seq[e+1]==12) && (seq[e]!=5) && (seq[e]!=12))
				pp[e] = true;
			else 
				pp[e] = false;
		pp[seq.length-1]=false;

		if (corpus.iSiteLib==null) {
		
		
/*		// Random picking
		for (int c=0; c<corpus.ungapped.length ; c++) {
			energies[c] = Math.random();
		}
*/		

/* Use this part of the code for the 'rmsChosen' scenario		
// Finding the fragment in the library
		int matchRes = -1;
		for (int c=0; c<corpus.ungapped.length ; c++) {
			boolean match = true;
			for (int e=0; (e<seq.length) && match ; e++) {
				if (seq[e] != corpus.resType[corpus.ungapped[c]+e])
					match = false;
			}
			if (match) {
				System.out.println("XXXXXXX: " + corpus.proteinNames[corpus.protInd[corpus.ungapped[c]]] + " " + corpus.resNum[corpus.ungapped[c]]);
				matchRes = corpus.ungapped[c];
				break;
			}
		}
*/
		
		
// **** This is for the loops paper **************************		
 		// Calculating the propensity threading energy. Putting NaN for fragments not complying with the Pro,Gly switching.
		double[] statBlos = new double[energies.length];
		double[] statProp = new double[energies.length];
 		for (int c=0; c<corpus.ungapped.length ; c++) {
 			boolean badGPswitch = false;
			energies[c] = 0.0;
			statBlos[c] = 0.0;
			statProp[c] = 0.0;
//			for (int e=0; e<seq.length ; e++) {
			for (int e=(fragStartInSeq-7); e<=(fragEndInSeq+7) ; e++) {
				if (pp[e]) {
					energies[c] += (W_PROP_ENERGY*corpus.propensityAll.preProVal(corpus.torsions[corpus.ungapped[c]+e][1], corpus.torsions[corpus.ungapped[c]+e][2]));
					//energies[c] += (W_PROP_ENERGY*corpus.prePro[corpus.ungapped[c]+e]);
					statProp[c] +=  corpus.propensityAll.preProVal(corpus.torsions[corpus.ungapped[c]+e][1], corpus.torsions[corpus.ungapped[c]+e][2]);
				}
				else {
					//if ((seq[e]!=5) && (seq[e]!=12))
					energies[c] += (W_PROP_ENERGY*corpus.propensityAll.propVal(seq[e],corpus.torsions[corpus.ungapped[c]+e][1], corpus.torsions[corpus.ungapped[c]+e][2]));
					//energies[c] += (W_PROP_ENERGY*corpus.energies[corpus.ungapped[c]+e][0][seq[e]]);
					statProp[c] += corpus.propensityAll.propVal(seq[e],corpus.torsions[corpus.ungapped[c]+e][1], corpus.torsions[corpus.ungapped[c]+e][2]);
					/* This is the original version, which allowed other residues to replace pro and gly if their angles were right
					if ((seq[e]==12) && (corpus.resType[corpus.ungapped[c]+e]!=12) && (Math.abs(corpus.torsions[corpus.ungapped[c]+e][1]+1.0472)>0.2))
						badGPswitch = true;
					if ((seq[e]!=5) && (corpus.resType[corpus.ungapped[c]+e]==5) && (corpus.torsions[corpus.ungapped[c]+e][1]>0.0))
						badGPswitch = true; 
						End of the original version. */
					if ((statusPG==1) || ((statusPG==2) && (e>=fragStartInSeq) && (e<=fragEndInSeq))) {
						if ((seq[e]==12) && (corpus.resType[corpus.ungapped[c]+e]!=12))
							badGPswitch = true;
						if ((seq[e]!=5) && (corpus.resType[corpus.ungapped[c]+e]==5))
							badGPswitch = true;
					}
				}
				if ((e>=fragStartInSeq-1) && (e<=fragEndInSeq+1)) {
					energies[c] += (W_BLOSUM62_ENERGY*Corpus.blosum62[seq[e]][corpus.resType[corpus.ungapped[c]+e]]);
					statBlos[c] -= Corpus.blosum62[seq[e]][corpus.resType[corpus.ungapped[c]+e]];
				}
			}
			if (badGPswitch && (statusPG<3))
				energies[c] = Double.NaN;
			
// ****** End for the loop paper

/* Use this part of the code for the 'rmsChosen' scenario
			if (matchRes == corpus.ungapped[c])
				energies[c] = Double.NaN;
			else
				energies[c] = corpus.calcRmsBetweenStruct(matchRes, corpus.ungapped[c], seq.length, 1, seq.length);
*/
 		}
 		
// 		Arrays.sort(statBlos);
// 		double meanVal = 0.0;
// 		double meanXX = 0.0;
// 		for (int c=0 ; c<Math.min(statBlos.length,1000) ; c++) {
// 			meanVal += statBlos[c];
// 			meanXX += statBlos[c]*statBlos[c];
// 		}
// 		meanVal = meanVal/Math.min(statBlos.length,1000);
// 		meanXX = meanXX/Math.min(statBlos.length,1000);
// 		System.out.print("33333333 " + meanVal + " " + Math.sqrt(meanXX-meanVal*meanVal) + " ");
// 		
// 		Arrays.sort(statProp);
// 		meanVal = 0.0;
// 		meanXX = 0.0;
// 		for (int c=0 ; c<Math.min(statBlos.length,1000) ; c++) {
// 			meanVal += statProp[c];
// 			meanXX += statProp[c]*statProp[c];
// 		}
// 		meanVal = meanVal/Math.min(statProp.length,1000);
// 		meanXX = meanXX/Math.min(statProp.length,1000);
// 		System.out.println(meanVal + " " + Math.sqrt(meanXX-meanVal*meanVal) + " ");

 		
// 		 Adding a random twist...
// 		for (int c=0; c<energies.length ; c++) {
// 			if (!Double.isNaN(energies[c])) {
// 				energies[c] += 2.0*Math.random();
// 			}
// 		}
		}
		else {
//			double beta = 3.5;
//			double Wblos = 0.75;
//			System.out.println("Setting fragment library by Isite library. Beta is: " + beta);
//			double[] tmpEnergiesIsite = new double[energies.length];
//			double[] tmpEnergiesBlos = new double[energies.length];
//			corpus.iSiteLib.setThreadingEnergies(seq,fragStartInSeq,fragEndInSeq,tmpEnergiesIsite);
//            //double counter=0.0;
//            //double sum=0.0;
//            //double sum2=0.0;
//	 		for (int c=0; c<corpus.ungapped.length ; c++) {
//	 			tmpEnergiesBlos[c] = 0.0;
//	 			for (int e=(fragStartInSeq-1); e<=(fragEndInSeq+1) ; e++) {
//					tmpEnergiesBlos[c] += (-1.0*Corpus.blosum62[seq[e]][corpus.resType[corpus.ungapped[c]+e]]);
//				}	
//	 			//if (tmpEnergies[c]<50) {
//	 			//counter++;
//	 			//  sum += tmpEnergies[c];
//	 			//   sum2 += tmpEnergies[c]*tmpEnergies[c];
//	 			//                      }	 			
//	 		}
//	 		for (int c=0; c<corpus.ungapped.length ; c++) {
//	 			energies[c] = MeshiProgram.randomNumberGenerator().nextDouble() * 
//	 			Math.exp(beta* (Wblos*(tmpEnergiesBlos[c]-5.7)/5.15 + (1-Wblos)*(tmpEnergiesIsite[c]+1.11)/1.72)  );
//	 		}
//            //System.out.println("4 4 4 " + counter + " " + sum + " " + sum2);


			// This is for a boltzman potential:
			// ----------------------------------			
			double Wblos = 0.4;
			System.out.println("Wblos is: " + Wblos);
//			double[] tmpEnergiesIsite = new double[energies.length];
			double[] tmpEnergiesBlos = new double[energies.length];
			double[] tmpEnergiesIden = new double[energies.length];
			double[] tmpEnergiesProp = new double[energies.length];
//			corpus.iSiteLib.setThreadingEnergies(seq,fragStartInSeq,fragEndInSeq,tmpEnergiesIsite);
			for (int c=0; c<corpus.ungapped.length ; c++) {
				tmpEnergiesBlos[c] = 0.0;
				for (int e=(fragStartInSeq-1); e<=(fragEndInSeq+1) ; e++) {
					tmpEnergiesBlos[c] += (-1.0*Corpus.blosum62[seq[e]][corpus.resType[corpus.ungapped[c]+e]]);
				}	
				tmpEnergiesIden[c] = -1.0;
				if (manner==-1) {
					for (int e=(fragStartInSeq+overlap); e<=fragEndInSeq ; e++) {
						if (seq[e]!=corpus.resType[corpus.ungapped[c]+e]) {
							tmpEnergiesIden[c] = 0.0;
						}
					}
//					int tmpCounter = 0;
//					int seqStartInProt = corpus.iSiteLib.findSeqStartInProt(seq);		
//					DecimalFormat fmt = new DecimalFormat("0");
//					if (tmpEnergiesIden[c]<0.0) {
//						tmpCounter++;
//						System.out.print(tmpCounter + ".  ");
//						for (int e=(fragStartInSeq-1); e<=(fragEndInSeq+1) ; e++) {
//							System.out.print(Residue.nameOneLetter(corpus.resType[corpus.ungapped[c]+e])+" ");
//						}
//						System.out.print("   ");
//						for (int e=(fragStartInSeq-1); e<=(fragEndInSeq+1) ; e++) {
//							System.out.print(Residue.nameOneLetter(seq[e])+" ");
//						}
//						System.out.print("   ");
//						for (int e=(fragStartInSeq+overlap); e<=fragEndInSeq ; e++) {
//							System.out.print(" " + fmt.format(180.0/Math.PI*(corpus.torsions[corpus.ungapped[c]+e][1]-corpus.iSiteLib.phiPsiInSeq(seqStartInProt+e)[0])));
//							System.out.print(" " + fmt.format(180.0/Math.PI*(corpus.torsions[corpus.ungapped[c]+e][2]-corpus.iSiteLib.phiPsiInSeq(seqStartInProt+e)[1])) + " ");
//						}
//						System.out.println("   ");
//					}
				}
				else {
					for (int e=fragStartInSeq; e<=(fragEndInSeq-overlap) ; e++) {
						if (seq[e]!=corpus.resType[corpus.ungapped[c]+e]) {
							tmpEnergiesIden[c] = 0.0;
						}
					}
				}
				tmpEnergiesProp[c] = 0.0;
				for (int e=(fragStartInSeq-7); e<=(fragEndInSeq+7) ; e++) {
					if (pp[e]) {
						tmpEnergiesProp[c] += (corpus.propensityAll.preProVal(corpus.torsions[corpus.ungapped[c]+e][1], corpus.torsions[corpus.ungapped[c]+e][2]));
					}
					else {
						tmpEnergiesProp[c] += (corpus.propensityAll.propVal(seq[e],corpus.torsions[corpus.ungapped[c]+e][1], corpus.torsions[corpus.ungapped[c]+e][2]));
					}
				}	 			
//				if (tmpEnergiesIsite[c]<IsiteLib.veryHighEnergy/2.0) {
//					energies[c] = MeshiProgram.randomNumberGenerator().nextDouble() * Math.exp(beta*tmpEnergiesIsite[c]);
//					energies[c] = Wblos*tmpEnergiesBlos[c] + (1-Wblos)*tmpEnergiesIsite[c];
//				}
//				else {
//					energies[c] = MeshiProgram.randomNumberGenerator().nextDouble() * Math.exp(beta*(10.0 + tmpEnergiesBlos[c]));
//					energies[c] = 1000.0 + tmpEnergiesBlos[c];
//					tmpEnergiesIsite[c] = 1000.0 + tmpEnergiesBlos[c];
//				}
				
				if (tmpEnergiesIden[c] > -1.0) {
					energies[c] = Double.NaN;
				}
				else {
					energies[c] = tmpEnergiesBlos[c];
				}
				
//				energies[c] = 10000*tmpEnergiesIden[c] + tmpEnergiesBlos[c];
				energies[c] = tmpEnergiesBlos[c];
			}

//			double veryLargeInc = Math.pow(2,20);
//			double[] tmpEnergiesRandom = new double[energies.length];
//			for (int c=0 ; c<energies.length ; c++) {
//			    energies[c] = 2*energies.length*veryLargeInc + MeshiProgram.randomNumberGenerator().nextDouble()*veryLargeInc;
//			    tmpEnergiesRandom[c] = MeshiProgram.randomNumberGenerator().nextDouble();
//			}
//	        	int[] sortedBlos = AbstractLoopBuilder.findTopMinArray(tmpEnergiesBlos, 10*desiredLibSize, Double.MAX_VALUE);
//		        int[] sortedRand = AbstractLoopBuilder.findTopMinArray(tmpEnergiesProp, 10*desiredLibSize, Double.MAX_VALUE);
//			int blosCounter = 0;
//			int randCounter = 0;
//			for (int c=0 ; c<10*desiredLibSize ; c++) {
//			    if ((c % 2) == 1) { // do other than blos
//				if (energies[sortedBlos[randCounter]] > (energies.length*veryLargeInc)) {
//				    energies[sortedRand[randCounter]] = c*veryLargeInc;
//				}
//				randCounter++;
//			    }
//			    else {
//				if (energies[sortedBlos[blosCounter]] > (energies.length*veryLargeInc)) {
//				    energies[sortedBlos[blosCounter]] = c*veryLargeInc;
//				}
//				blosCounter++;
//			    }
//			}







			
			// This is for a combined potential:
			// ----------------------------------			
//			double Wblos = 0.99;
//			System.out.println("Setting fragment library by Isite library.");
//			double[] tmpEnergiesIsite = new double[energies.length];
//			double[] tmpEnergiesBlos = new double[energies.length];
//			corpus.iSiteLib.setThreadingEnergies(seq,fragStartInSeq,fragEndInSeq,tmpEnergiesIsite);
//			for (int c=0; c<corpus.ungapped.length ; c++) {
//	 			tmpEnergiesBlos[c] = 0.0;
//	 			for (int e=(fragStartInSeq-1); e<=(fragEndInSeq+1) ; e++) {
//					tmpEnergiesBlos[c] += (-1.0*Corpus.blosum62[seq[e]][corpus.resType[corpus.ungapped[c]+e]]);
//				}	
//	 			if (tmpEnergiesIsite[c]<IsiteLib.veryHighEnergy/2.0) {
//	 				energies[c] = Wblos*tmpEnergiesBlos[c] + (1-Wblos)*tmpEnergiesIsite[c];
//	 			}
//	 			else {
//	 				energies[c] = tmpEnergiesBlos[c];
//	 			}
//			}
			
			

// This is for an RMS-like potential:
// ----------------------------------			
//			int seqStartInProt = corpus.iSiteLib.findSeqStartInProt(seq);
//			System.out.println("Setting fragment library by Isite library. Frag start in residue: " + seqStartInProt);
//	 		for (int c=0; c<corpus.ungapped.length ; c++) {
//				double maxDevBB = 0.0;
//	 			for (int res=fragStartInSeq; res<=fragEndInSeq ; res++) {
//					double tmpDev = IsiteLib.torsionDiff(corpus.torsions[corpus.ungapped[c]+res][1],
//							corpus.iSiteLib.phiPsiInSeq(seqStartInProt+res)[0]);
//					if (tmpDev>maxDevBB) {
//						maxDevBB = tmpDev;
//					}
//					tmpDev = IsiteLib.torsionDiff(corpus.torsions[corpus.ungapped[c]+res][2],
//							corpus.iSiteLib.phiPsiInSeq(seqStartInProt+res)[1]);
//					if (tmpDev>maxDevBB) {
//						maxDevBB = tmpDev;
//					}
//				}
//	 			energies[c] = maxDevBB;
//			}
		
		
		}
 		return energies;
	}

	
	protected void buildLib(boolean filterStalks) {
		
		// Initialization:
		int[] tmpIndices = new int[desiredLibSize];
		double[] tmpEnergies = new double[desiredLibSize];
		for (int c=0; c<desiredLibSize ; c++) {
			tmpIndices[c] = -1;
			tmpEnergies[c] = Double.MAX_VALUE;
		}

		// Start iterating for the lib
		int currentLibSize=0;
		for ( ; currentLibSize<desiredLibSize ; currentLibSize++) {
			int lowInd = findLowestEnergy(tmpIndices,threadingEnergies,filterStalks,fragStartInSeq,currentLibSize);
			if (lowInd==-1) {
				issueSmallLibWarning(currentLibSize);
				libSize=currentLibSize;
				break;
			}
			tmpIndices[currentLibSize] = corpus.ungapped[lowInd] + fragStartInSeq;
			tmpEnergies[currentLibSize] = threadingEnergies[lowInd];
			similarToLib[lowInd] = true;
//			if (currentLibSize == (3*desiredLibSize/4)) {
//				for (int fragCounter=0; fragCounter<threadingEnergies.length ; fragCounter++) {
//					threadingEnergies[fragCounter] = MeshiProgram.randomNumberGenerator().nextDouble() * Math.exp(0.1*threadingEnergies[fragCounter]);
//				}
//				System.out.print("Went to boltzman rep.\n");
//			}
		}
		libSize = currentLibSize;

		// Updating the lib fields according to its final size
		libOrig = new int[libSize];
		libEnergy =  new double[libSize];
		for (int c=0 ; c<libSize ; c++) {
			libOrig[c] = tmpIndices[c];
			libEnergy[c] = tmpEnergies[c];
		}
		
		
		// Just ploting the max deviation
		DecimalFormat fmt = new DecimalFormat("0.#");
		int seqStartInProt = corpus.iSiteLib.findSeqStartInProt(seq);
		for (int c=0 ; c<libSize ; c++) {
			System.out.print("3 3 3 " + c);
			if (manner==-1) {
				for (int res=0; res<fragL ; res++) {
					System.out.print(" " + fmt.format(180.0/Math.PI*IsiteLib.torsionDiff(corpus.torsions[libOrig[c]+res][1],corpus.iSiteLib.phiPsiInSeq(seqStartInProt+res)[0])));
					System.out.print(" " + fmt.format(180.0/Math.PI*IsiteLib.torsionDiff(corpus.torsions[libOrig[c]+res][2],corpus.iSiteLib.phiPsiInSeq(seqStartInProt+res)[1])) + " ");
				}
			}
			else if (manner==1) {
				for (int res=(fragL-1); res>-1 ; res--) {
					System.out.print(" " + fmt.format(180.0/Math.PI*IsiteLib.torsionDiff(corpus.torsions[libOrig[c]+res][1],corpus.iSiteLib.phiPsiInSeq(seqStartInProt+res)[0])));
					System.out.print(" " + fmt.format(180.0/Math.PI*IsiteLib.torsionDiff(corpus.torsions[libOrig[c]+res][2],corpus.iSiteLib.phiPsiInSeq(seqStartInProt+res)[1])) + " ");
				}
			}
			System.out.println();
		}
		
	}  // Of build lib


	public String printLib() {
		//if (calcRmsBetweenStruct(libOrig[0],corpus.ungapped[187]) < 1.0)
		//	return "Library discarded. Too close to a canonical helix.\n";
		String result = "";

		for (int i=0; i<libSize ; i++) {
			result += i + " " + libEnergy[i] + " " + 
			corpus.proteinNames[corpus.protInd[libOrig[i]]] + " " + corpus.resNum[libOrig[i]] + "\n";
		}

		return result;
	}

	public String printLib(int i) {
		//if (calcRmsBetweenStruct(libOrig[0],corpus.ungapped[187]) < 1.0)
		//	return "Library discarded. Too close to a canonical helix.\n";
		String result = "";

		result += i + " " + libEnergy[i] + " " + 
		corpus.proteinNames[corpus.protInd[libOrig[i]]] + " " + corpus.resNum[libOrig[i]] + "\n";

		return result;
	}


	/** This method checks if the fragment 'ind' (given in the ungapped indexing) is similar to a fragment already in the lib.
 The array similarInLib is updated if a similarity is present.  
	 **/
	private boolean checkSimilarity(int ind, int[] tmpIndices, int fragStartInSeq, boolean filterStalks) {
		if (similarToLib[ind])
			return true;
		for (int i=0; tmpIndices[i]>-1 ; i++) 
			if (!filterStalks) {
				if (corpus.calcRmsBetweenStruct(tmpIndices[i],corpus.ungapped[ind]+fragStartInSeq,fragL,1,fragL) < rmsCutOff) {
					similarToLib[ind] = true;
					return true;
				}
			}
			else {
				if (corpus.calcRmsBetweenStruct(tmpIndices[i],corpus.ungapped[ind]+fragStartInSeq,fragL,manner,overlap) < rmsCutOff) {
					similarToLib[ind] = true;
					return true;
				}
			}
		return false;
	}

	/** This method finds the lowest energy fragment that is not too similar to the rest of the lib **/
	private int findLowestEnergy(int[] tmpIndices, double[] energies, boolean filterStalks, int fragStartInSeq,int currentLibSize) {
		double lowestEnergy = Double.MAX_VALUE;
		int lowestEnergyInd = -1;
		do {
			lowestEnergy = Double.MAX_VALUE;
			lowestEnergyInd = -1;
			for (int i=0; i<corpus.ungapped.length ; i++) {
				if ( (energies[i]<lowestEnergy) && !similarToLib[i] ) {
					if (!filterStalks /* || (currentLibSize < 0.5*desiredLibSize) */ || isFragCompatibleForConstruction(corpus.ungapped[i] + fragStartInSeq , overlapSimilarityTH)) {
						if ((prevLib==null) || 
								isFragCompatibleToPrevLib(corpus.ungapped[i] + fragStartInSeq, 
										prevLib.libSize*currentLibSize/desiredLibSize, 
										overlap, manner, overlapSimilarityTH)) { // I added this on 4/26/2010
							lowestEnergy = energies[i];
							lowestEnergyInd = i;
						}
					}
					else {
						similarToLib[i] = true;
					}
				}
			}
		} while ((lowestEnergyInd>-1) && checkSimilarity(lowestEnergyInd,tmpIndices,fragStartInSeq,filterStalks));
		return lowestEnergyInd;
	}


	/** 
This method check whether the fragment 'indInLib' in the lib can overlap the stalks below a certain
torsion-diff cut off. The fragment correspond to a protein fragment starting in residue 'residueFragStart'.
The protein torsions themself is taken from the Isite lib.
'overlap' is the number of residues that are used for assessing torsion-diff. 'manner' means:
 -1  - overlap of the N-terminus
 0 - Not implemented and will give a runtime exception. (two overlap regions in both termini (each of length overlap))
 1  - overlap of the C-terminus
 
 **/
 
	public boolean isFragCompatible(int indInLib, int residueFragStart, int overlap, int manner, double overlapTH) {
		return isFragCompatibleCorpusIndex(libOrig[indInLib], residueFragStart, overlap, manner, overlapTH);
	} 

	protected boolean isFragCompatibleCorpusIndex(int indInCorpus, int residueFragStart, int overlap, int manner, double overlapTH) {
		if (manner==-1) {
			for (int c=0 ; c<overlap ; c++) {
				if(IsiteLib.torsionDiff(corpus.torsions[indInCorpus+c][1],
						corpus.iSiteLib.phiPsiInSeq(residueFragStart+c)[0])>overlapTH) {
					return false;
				}
				if(IsiteLib.torsionDiff(corpus.torsions[indInCorpus+c][2],
						corpus.iSiteLib.phiPsiInSeq(residueFragStart+c)[1])>overlapTH) {
					return false;
				}
			}
			return true;
		}
		if (manner==1) {
			for (int c=fragL-overlap ; c<fragL ; c++) {
				if(IsiteLib.torsionDiff(corpus.torsions[indInCorpus+c][1],
						corpus.iSiteLib.phiPsiInSeq(residueFragStart+c)[0])>overlapTH) {
					return false;
				}
				if(IsiteLib.torsionDiff(corpus.torsions[indInCorpus+c][2],
						corpus.iSiteLib.phiPsiInSeq(residueFragStart+c)[1])>overlapTH) {
					return false;
				}
			}
			return true;
		}
		throw new RuntimeException("Unrecognized manner: " + manner);
	} 
	
	/**
	 *The 'overlap', 'manner' and 'residueFragStart' are from the class's fields. 
	 **/
	public boolean isFragCompatible(int indInLib, double overlapTH) {
		return isFragCompatibleCorpusIndex(libOrig[indInLib], residueFragStart, overlap, manner,  overlapTH);
	} 	
	
	/**
	 *This method is used only during the construction when 'libOrig' is not yet ready. 
	 **/
	protected boolean isFragCompatibleForConstruction(int indInCorpus, double overlapTH) {
		if (libOrig!=null)
			throw new RuntimeException("\nERROR to the programer: you should call this method only if libOrig is still unavailable (i.e. ==null) because the constructor is running. \n");
		return isFragCompatibleCorpusIndex(indInCorpus, residueFragStart, overlap, manner,  overlapTH);
	} 	
	

	protected boolean isFragCompatibleToPrevLib(int indInCorpus, int indInOtherLib, int overlap, int manner, double overlapTH) {
		return true;
		
		
		
//		int indInCorpusPrevLib = prevLib.libOrig(indInOtherLib);
//		int fragLprevLib = prevLib.fragL();
//		if (manner==-1) {
//			for (int c=0 ; c<overlap ; c++) {
//				if(IsiteLib.torsionDiff(corpus.torsions[indInCorpus+c][1],
//						corpus.torsions[indInCorpusPrevLib+fragLprevLib-overlap+c][1])>overlapTH) {
//					return false;
//				}
//				if(IsiteLib.torsionDiff(corpus.torsions[indInCorpus+c][2],
//						corpus.torsions[indInCorpusPrevLib+fragLprevLib-overlap+c][2])>overlapTH) {
//					return false;
//				}
//			}
//			return true;
//		}
//		if (manner==1) {
//			for (int c=fragL-overlap ; c<fragL ; c++) {
//				if(IsiteLib.torsionDiff(corpus.torsions[indInCorpus+c][1],
//						corpus.torsions[indInCorpusPrevLib+(c-(fragL-overlap))][1])>overlapTH) {
//					return false;
//				}
//				if(IsiteLib.torsionDiff(corpus.torsions[indInCorpus+c][2],
//						corpus.torsions[indInCorpusPrevLib+(c-(fragL-overlap))][2])>overlapTH) {
//					return false;
//				}
//			}
//			return true;
//		}
//		throw new RuntimeException("Unrecognized manner: " + manner);

	} 
	
	
	
	/**
	 * Will insert the fragment with that index into the protein.
	 * The 'prot', 'overlap', 'manner' and 'residueFragStart' are from the class's fields. 
	 **/ 
	public void insertFragToProt(DunbrackLib rotLib, int indInLib) {
			insertFragToProt(rotLib, indInLib, prot, residueFragStart, overlap, manner);
		}
	
	/**
	 * Will insert the fragment with that index into the protein.
	 * 'prot' neednot be the protein from which the sequence for the library was taken.
	 **/ 
	public void insertFragToProt(DunbrackLib rotLib, int indInLib, Protein prot, int residueFragStart, int overlap, int manner) {
		if (manner==-1) {  // N term
			Residue prevRes = prot.residue(residueFragStart+overlap-1);
			for (int c=overlap ; c<fragL ; c++) {
				Residue buildRes = prot.residue(residueFragStart+c);
				ResidueBuilder.buildBackboneFromNterm(prevRes,
						buildRes, 
						corpus.torsions[libOrig[indInLib]+c-1][2], 
						corpus.torsions[libOrig[indInLib]+c][0],
						corpus.torsions[libOrig[indInLib]+c][1],
						corpus.torsions[libOrig[indInLib]+c][2]);
				double[] rot = null;
				if ((seq[c]!=0)&&(seq[c]!=5))
					rot = rotLib.getRotamer(seq[c], corpus.torsions[libOrig[indInLib]+c][1],
							corpus.torsions[libOrig[indInLib]+c][2] , 0);
				ResidueBuilder.build(buildRes, buildRes.type, rot); 
				corpus.iSiteLib.phiPsiInSeq(residueFragStart+c)[0] = corpus.torsions[libOrig[indInLib]+c][1];
				corpus.iSiteLib.phiPsiInSeq(residueFragStart+c)[1] = corpus.torsions[libOrig[indInLib]+c][2];
				prevRes = buildRes;
			}
		}
		else if (manner==1) {  // C term
			Residue postRes = prot.residue(residueFragStart+fragL-overlap);
			for (int c=fragL-overlap-1 ; c>-1 ; c--) {
				Residue buildRes = prot.residue(residueFragStart+c);
				ResidueBuilder.buildBackboneFromCterm(postRes,
						buildRes, 
						corpus.torsions[libOrig[indInLib]+c][1],
						corpus.torsions[libOrig[indInLib]+c][2],
						corpus.torsions[libOrig[indInLib]+c+1][0],
						corpus.torsions[libOrig[indInLib]+c+1][1]);
				double[] rot = null;
				if ((seq[c]!=0)&&(seq[c]!=5))
					rot = rotLib.getRotamer(seq[c], corpus.torsions[libOrig[indInLib]+c][1],
							corpus.torsions[libOrig[indInLib]+c][2] , 0);
				ResidueBuilder.build(buildRes, buildRes.type, rot); 
				corpus.iSiteLib.phiPsiInSeq(residueFragStart+c)[0] = corpus.torsions[libOrig[indInLib]+c][1];
				corpus.iSiteLib.phiPsiInSeq(residueFragStart+c)[1] = corpus.torsions[libOrig[indInLib]+c][2];
				postRes = buildRes;
			}			
		}
		else {  // both sides
			throw new RuntimeException("\n\nNot yet implemented - fragment insertion restrained on both sides. Im not even sure if it can be defined in torsion coordinates.\n\n");
		}
	}

	/**
	 *  A warning issued if the lib failed to reach the desired size.
	 */  
	protected void issueSmallLibWarning(int currentLibSize) {
		System.out.println("WARNINIG!!! The corpus is too small to create such a large library. So Far created " + currentLibSize);
	}

	/**
	 * Any entry that is "true" in the similarToLib array, will cause the corresponding fragment from the corpus
	 * NOT to be considered in the lib. This means that manipulating the entries in this array is a great way to
	 * direct which fragments are optional from the corpus.
	 */
	protected void initializeSimilarToLib() {
		similarToLib = new boolean[corpus.ungapped.length];
		for (int c=0; c<corpus.ungapped.length ; c++)
			similarToLib[c] = false;	
	}
	

//	Getters...
	protected Corpus corpus() {return corpus;}   // For the inheriting classes
	protected Protein prot() {return prot;}   // For the inheriting classes
	protected int fragL() {return fragL;}   // For the inheriting classes
	public int manner() {return manner;}   
	public int overlap() {return overlap;}   
	public double libEnergy(int ind) { return libEnergy[ind];}
	public int libSize() { return libSize;}
	public int libOrig(int ind) { return libOrig[ind];}
	public int trueFragStart() {
		if (prot==null)
			throw new RuntimeException("\nThis methods can be run only for instances that were initiated with constructor No. 1\n");
		else
			return trueFragStart; 
	}
	public int trueFragEnd() {		
		if (prot==null)
			throw new RuntimeException("\nThis methods can be run only for instances that were initiated with constructor No. 1\n");
		else
			return trueFragEnd; 
	}
	public int residueFragStart() { 
		if (prot==null)
			throw new RuntimeException("\nThis methods can be run only for instances that were initiated with constructor No. 1\n");
		else
			return residueFragStart; 
	}


} // Of LocalFragLib
