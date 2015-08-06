package meshi.applications.localFragLibs;

import meshi.applications.corpus.Corpus;
import meshi.molecularElements.Protein;
import meshi.util.Isite.IsiteLib;

/**
 * 
 * @author Nir
 *
 * The purpose of this class is to create a library with only one segment that should be the
 * closest fragment in the corpus to a certain part of a given protein. 
 * 
 */


public class LocalSingleFragmentLib extends LocalFragmentLib {

	public LocalSingleFragmentLib(Corpus corpus, int[] seq, int libSize,
			double rmsCutoff, int fragStartInSeq, int fragEndInSeq,
			Protein prot, int residueFragStart, int overlap, int manner,
			double overlapRMSsimilarity, int statusPG) {
		super(corpus, seq, libSize, rmsCutoff, fragStartInSeq, fragEndInSeq,
				prot, residueFragStart, overlap, manner, overlapRMSsimilarity, statusPG);
	}

	/**
	 * As this fragment does not need to overlap with the protein, we want to pick it always.
	 **/
	protected boolean isFragCompatibleForConstruction(int indInCorpus, double overlapRMSsimilarity) {
		return true;
	} 	

	
	
	protected void initializeSimilarToLib() {
		similarToLib = new boolean[corpus().ungapped.length];
		for (int c=0; c<corpus().ungapped.length ; c++)
			similarToLib[c] = true;
		similarToLib[findBestFittingFrag()] = false;
	}
	
	private int findBestFittingFrag() {
		double bestFitSoFar = 1e100;
		int bestFrag = -1;
		for (int c=corpus().ungapped.length-1; c>-1 ; c--) {  // Going from last, because the model is there.
			//double tmpFit = 0; // superimposeFrag(corpus().ungapped[c], prot(), residueFragStart(), fragL(), manner());
			double tmpFit = 0.0;
			for (int res=0; res<fragL() ; res++) {
				double tmpDev = IsiteLib.torsionDiff(corpus().torsions[corpus().ungapped[c]+res][1],
						corpus().iSiteLib.phiPsiInSeq(residueFragStart()+res)[0]);
				if (tmpDev>tmpFit) {
					tmpFit = tmpDev;
				}
				tmpDev = IsiteLib.torsionDiff(corpus().torsions[corpus().ungapped[c]+res][2],
						corpus().iSiteLib.phiPsiInSeq(residueFragStart()+res)[1]);
				if (tmpDev>tmpFit) {
					tmpFit = tmpDev;
				}			
			}			
			if (tmpFit<bestFitSoFar) {
				bestFitSoFar = tmpFit;
				bestFrag = c;
			}
			if (bestFitSoFar<0.001) // If you found a very good frag, don't waste time.
				return bestFrag;
		}
		return bestFrag;
	}

}
