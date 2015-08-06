package meshi.applications.localFragLibs;

import meshi.applications.corpus.Corpus;
import meshi.molecularElements.Protein;

/**
 * 
 * @author Nir
 *
 * This class is intended to build a fragment lib that on one side interface with a peptide chain, and on the other side interface 
 * a (single) fragment from the LocalSingleFragmentLib. This means that the fragments in this lib must be overlapping consistent with
 * that single fragment. For this reason a "bridging" lib (i.e. manner = 0) is not allowed. The way this overlapping consistency is  
 * achieved is by filtering the corpus for this property with the "initializeSimilarToLib" method. Because we don't want a lib to be 
 * created during the constructor's "super" call, the "initializeSimilarToLib" method gives 0 possible fragments until the local field:
 * "singleLib" has been set.
 *  
 */
public class LocalFragmentLibInterfacingToSingle extends LocalFragmentLib {

	private LocalSingleFragmentLib singleLib=null;
	
	public LocalFragmentLibInterfacingToSingle(Corpus corpus, int[] seq,
			int libSize, double rmsCutoff, int fragStartInSeq,
			int fragEndInSeq, Protein prot, int residueFragStart, int overlap,
			int manner, double overlapSimilarityTH, LocalSingleFragmentLib singleLib, int statusPG) {
		super(corpus, seq, libSize, rmsCutoff, fragStartInSeq, fragEndInSeq,
				prot, residueFragStart, overlap, manner, overlapSimilarityTH, statusPG);
		if (manner==0)
			throw new RuntimeException("This class could only work with manner set to '1' or '-1'. See documentation why this is so.");
		this.singleLib = singleLib;
		buildLib(true);
	}

	public LocalFragmentLibInterfacingToSingle(Corpus corpus, int[] seq,
			int libSize, double rmsCutoff, int fragStartInSeq, int fragEndInSeq, 
			int manner, double overlapSimilarityTH, int overlap, LocalSingleFragmentLib singleLib, int statusPG) {
		super(corpus, seq, libSize, rmsCutoff, fragStartInSeq, fragEndInSeq, manner, overlap, statusPG, 0.00000, null);
		if (true)
			throw new RuntimeException("On 4/26/2010 I change the constructor of LocalFragLib to get overlap_similarity and the previous_lib. As a result I needed to add the 0.0 and null in the super two lines above, but I guess it will not work now.");
		this.singleLib = singleLib;
		this.overlapSimilarityTH = overlapSimilarityTH;
		buildLib(false);
	}

	/**
	 *  Only on the second build of the library we want to issue a warning.
	 */  
	protected void issueSmallLibWarning(int currentLibSize) {
		if (singleLib != null)
			System.out.println("WARNINIG!!! The corpus is too small to create such a large library. So Far created " + currentLibSize);
	}

	/**
	 * We discard the first build of the library in the constructor's "super" call by putting all values to "true".
	 */
	protected void initializeSimilarToLib() {
		if (singleLib == null) {// First dummy built of the library, all values are true and no fragment will be taken from the corpus.
			similarToLib = new boolean[corpus().ungapped.length];
			for (int c=0; c<corpus().ungapped.length ; c++)
				similarToLib[c] = true;	
		}
		else { // The second important pass.
			similarToLib = new boolean[corpus().ungapped.length];
			if (manner()==-1) { // This lib is connecting frags to a Cterm of the chain, and so the single lib is interfacing to the Cterm of the frags.
				for (int c=0; c<corpus().ungapped.length ; c++) 
					if (corpus().calcRmsBetweenStruct(corpus().ungapped[c]+fragL()-overlap(), singleLib.libOrig(0), overlap(), 1, overlap())<overlapSimilarityTH)
						similarToLib[c] = false;
					else
						similarToLib[c] = true;
			}
			else { // This lib is connecting frags to a Nterm of the chain, and so the single lib is interfacing to the Nterm of the frags.
				for (int c=0; c<corpus().ungapped.length ; c++) 
					if (corpus().calcRmsBetweenStruct(corpus().ungapped[c], singleLib.libOrig(0)+singleLib.fragL()-overlap(), overlap(), 1, overlap())<overlapSimilarityTH)
						similarToLib[c] = false;
					else
						similarToLib[c] = true;				
			}
		}
	}

	
	/**
	 * This override allows the method to be called also on the second time. 
	 **/
	protected boolean isFragCompatibleForConstruction(int indInCorpus, double overlapTH) {
		return isFragCompatibleCorpusIndex(indInCorpus, residueFragStart(), overlap(), manner(),  overlapTH);
	} 	

}
