package meshi.applications.loopBuilding;

import java.util.Vector;

import meshi.applications.corpus.Corpus;
import meshi.applications.localFragLibs.LocalFragmentLib;
import meshi.applications.localFragLibs.LocalFragmentLibInterfacingToSingle;
import meshi.applications.localFragLibs.LocalSingleFragmentLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;

/** 
 * This class is building the local library array according the user input. 
 *  
 * @author Nir
 *
 */

public abstract class AbstractLoopBuilderUserDefinedFragments extends AbstractLoopBuilder {

	protected final double EV_CUTOFF_FACTOR=1.5; // The building tree is pruned if the EV is above this value * loop length 
	private double tetherToTemplate = -1; // If this is set to anything larger than 0, then the fragment elongation will not proceed if the inserted fragment CA's are more than this values from theirs initial-structure counterparts. 
	private double tetherToTemplateSquared = -1; // The square of the above.

	protected Vector<int[]> fragsDescription; 
	protected double rmsCutOff = -1;
	protected double EV_CUTOFF = 0.0;
	protected long EVaccepted = 0;
	protected long EVrejected = 0;
	protected long TETHERaccepted = 0;
	protected long TETHERrejected = 0;
	private Atom[] templateCAs = null;
	private Atom[] modelCAs = null;
	
	/**
	 * 
	 * New parameters not in the 'super':
	 * @param rmsCutOff - that should be the RMS difference between the fragments in the library
	 * @param fragsDescription - This is a vector of 4-tupples. In each 4 ints the data is: [start,end,type,lib size]
	 * The numbering of the type are:
	 *  -1,1 - regular N- C- term frags, connecting to protein.
	 *  -2,2 - like above but not connecting to a known structure.
	 *  -3,3 - like above but fixed.
	 *   <-999, >999 , like above but interfacing a fixed frags on the other side. The frag index is found by (X-1000) or (-1000-X) depending on the sign. 
	 * 
	 */
	public AbstractLoopBuilderUserDefinedFragments(CommandList commands,
			String writePath, Corpus corpus, Protein prot, Protein ref,
			int resStart, int resEnd,
			double similarityMatchCO, double rmsCutOff, Vector<int[]> fragsDescription) {
		super(commands, writePath, corpus, prot, ref, resStart, resEnd,
				similarityMatchCO);
		this.fragsDescription = fragsDescription;
		this.rmsCutOff = rmsCutOff;
		EV_CUTOFF = EV_CUTOFF_FACTOR*(resEnd-resStart+1);
		System.out.println("EV_CUTOFF: " + EV_CUTOFF);	
		setTemplateTetherStructure();
	}
	
	
	// Build the local frag libraries.
	protected void buildLibraries() {
		// The auxillary arrays
		libs = new LocalFragmentLib[fragsDescription.size()];
		libStarts = new int[fragsDescription.size()];    // The indices in the loop (first missing res is 0) where the fragsDescription (not overlaps) of each lib strat.
		libEnds = new int[fragsDescription.size()];      // The same as above, but where they end
		libManners = new int[fragsDescription.size()];  // The manner of each lib. (-1) - overlapping in the N terminus. (1) - overlapping in the C terminus.
		startResForClosingPotential = new int[fragsDescription.size()]; // For each lib, this and the next array can be used as the parameters for the term 'temporaryDistEnergy'
		endResForClosingPotential = new int[fragsDescription.size()];
		fragRank = new int[fragsDescription.size()];
		actuallyTaken = new int[fragsDescription.size()];
		// First building the fixed libraries.
		for (int c=0 ; c<fragsDescription.size() ; c++) {
			if (fragsDescription.get(c)[2]==-3) {
				libs[c] = new LocalSingleFragmentLib(corpus, 
						Protein.getSeqOfProt(prot, fragsDescription.get(c)[0], fragsDescription.get(c)[1]), 
						1, 99.9, 0, fragsDescription.get(c)[1]-fragsDescription.get(c)[0], prot, fragsDescription.get(c)[0],
						overlapStalk, -1, similarityMatchCO, 3 /*No PG check is needed here*/);
				libStarts[c] = fragsDescription.get(c)[0]-resStart+overlapStalk;
				libEnds[c] = fragsDescription.get(c)[1]-resStart;
				libManners[c] = -1;
			}
			if (fragsDescription.get(c)[2]==3) {
				libs[c] = new LocalSingleFragmentLib(corpus, 
						Protein.getSeqOfProt(prot, fragsDescription.get(c)[0], fragsDescription.get(c)[1]), 
						1, 99.9, 0, fragsDescription.get(c)[1]-fragsDescription.get(c)[0], prot, fragsDescription.get(c)[0],
						overlapStalk, 1, similarityMatchCO, 3 /*No PG check is needed here*/);
				libStarts[c] = fragsDescription.get(c)[0]-resStart;
				libEnds[c] = fragsDescription.get(c)[1]-overlapStalk-resStart;
				libManners[c] = 1;
			}
		}
		// Next building the library as they appear in the Vector
		for (int c=0 ; c<fragsDescription.size() ; c++) {
			if (fragsDescription.get(c)[2]<-999) {
				int FRAGSIZE = fragsDescription.get(c)[1]-fragsDescription.get(c)[0]+1+overlapStalk;
				libs[c] = new LocalFragmentLibInterfacingToSingle(corpus , Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading,fragsDescription.get(c)[0]+(FRAGSIZE-1)) , 
							fragsDescription.get(c)[3],  rmsCutOff, (extensionThreading-overlapStalk) , (extensionThreading+(FRAGSIZE-1)) ,  
							prot,(fragsDescription.get(c)[0]-overlapStalk),overlapStalk,-1,similarityMatchCO, (LocalSingleFragmentLib) libs[-1000-fragsDescription.get(c)[2]],3 /*No PG check is needed here*/);
				libStarts[c] = fragsDescription.get(c)[0]-resStart;
				libEnds[c] = libStarts[c] + (FRAGSIZE-1);
				libManners[c] = -1;
				if (true) { // To print
					System.out.println("FRAGSIZE " + FRAGSIZE);
					for (int d=0; d<Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading,fragsDescription.get(c)[0]+(FRAGSIZE-1)).length ; d++)
						System.out.print(Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading,fragsDescription.get(c)[0]+(FRAGSIZE-1))[d]+ " ");
					System.out.println();
					System.out.println(libStarts[c] + "   SSS   " + libEnds[c]);
					System.out.println("variable Lib Size is:" + libs[c].libSize());
				}
			}
			if (fragsDescription.get(c)[2]>999) {
				int FRAGSIZE = fragsDescription.get(c)[1]-fragsDescription.get(c)[0]+1+overlapStalk;
				libs[c] = new LocalFragmentLibInterfacingToSingle(corpus , Protein.getSeqOfProt(prot, fragsDescription.get(c)[1]-(FRAGSIZE-1),fragsDescription.get(c)[1]+extensionThreading), 
							fragsDescription.get(c)[3],  rmsCutOff, 0 , (FRAGSIZE-1)+overlapStalk ,  
							prot,(fragsDescription.get(c)[1]-(FRAGSIZE-1)),overlapStalk, 1 ,similarityMatchCO, (LocalSingleFragmentLib) libs[fragsDescription.get(c)[2]-1000],3 /*No PG check is needed here*/);
				libStarts[c] = fragsDescription.get(c)[0] - overlapStalk - resStart;
				libEnds[c] = libStarts[c] + (FRAGSIZE-1);
				libManners[c] = 1;
				if (true) { // To print
					System.out.println("FRAGSIZE " + FRAGSIZE);
					for (int d=0; d<Protein.getSeqOfProt(prot, fragsDescription.get(c)[1]-(FRAGSIZE-1),fragsDescription.get(c)[1]+extensionThreading).length ; d++)
						System.out.print(Protein.getSeqOfProt(prot, fragsDescription.get(c)[1]-(FRAGSIZE-1),fragsDescription.get(c)[1]+extensionThreading)[d]+ " ");
					System.out.println();
					System.out.println(libStarts[c] + "   SSS   " + libEnds[c]);
					System.out.println("variable Lib Size is:" + libs[c].libSize());
				}
			}
			if (fragsDescription.get(c)[2]==-1) {
				int FRAGSIZE = fragsDescription.get(c)[1]-fragsDescription.get(c)[0]+1;
				libs[c] = new LocalFragmentLib(corpus , Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading,fragsDescription.get(c)[0]+(FRAGSIZE-1)+extensionThreading) , 
							fragsDescription.get(c)[3],  rmsCutOff, (extensionThreading-overlapStalk) , (extensionThreading+(FRAGSIZE-1)) ,  
							prot,(fragsDescription.get(c)[0]-overlapStalk),overlapStalk,-1,similarityMatchCO,3 /*No PG check is needed here*/);
				libStarts[c] = fragsDescription.get(c)[0] - resStart;
				libEnds[c] = libStarts[c] + (FRAGSIZE-1);
				libManners[c] = -1;
				if (true) { // To print
					System.out.println("FRAGSIZE " + FRAGSIZE);
					for (int d=0; d<Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading,fragsDescription.get(c)[0]+(FRAGSIZE-1)+extensionThreading).length ; d++)
						System.out.print(Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading,fragsDescription.get(c)[0]+(FRAGSIZE-1)+extensionThreading)[d]+ " ");
					System.out.println();
					System.out.println(libStarts[c] + "   SSS   " + libEnds[c]);
					System.out.println("variable Lib Size is:" + libs[c].libSize());
				}
			}
			if (fragsDescription.get(c)[2]==1) {
				int FRAGSIZE = fragsDescription.get(c)[1]-fragsDescription.get(c)[0]+1;
				libs[c] = new LocalFragmentLib(corpus , Protein.getSeqOfProt(prot, fragsDescription.get(c)[1]-(FRAGSIZE-1)-extensionThreading,fragsDescription.get(c)[1]+extensionThreading) , 
							fragsDescription.get(c)[3],  rmsCutOff, extensionThreading , (extensionThreading+(FRAGSIZE-1)+overlapStalk) ,  
							prot,(fragsDescription.get(c)[1]-(FRAGSIZE-1)),overlapStalk,1,similarityMatchCO,3 /*No PG check is needed here*/);
				libStarts[c] = fragsDescription.get(c)[0] - resStart;
				libEnds[c] = libStarts[c] + (FRAGSIZE-1);
				libManners[c] = 1;
				if (true) { // To print
					System.out.println("FRAGSIZE " + FRAGSIZE);
					for (int d=0; d<Protein.getSeqOfProt(prot, fragsDescription.get(c)[1]-(FRAGSIZE-1)-extensionThreading,fragsDescription.get(c)[1]+extensionThreading).length ; d++)
						System.out.print(Protein.getSeqOfProt(prot, fragsDescription.get(c)[1]-(FRAGSIZE-1)-extensionThreading,fragsDescription.get(c)[1]+extensionThreading)[d]+ " ");
					System.out.println();
					System.out.println(libStarts[c] + "   SSS   " + libEnds[c]);
					System.out.println("variable Lib Size is:" + libs[c].libSize());
				}
			}
			if (fragsDescription.get(c)[2]==-2) {
				int FRAGSIZE = fragsDescription.get(c)[1]-fragsDescription.get(c)[0]+1;
				libs[c] = new LocalFragmentLib(corpus , Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading, fragsDescription.get(c)[0]+(FRAGSIZE-1)+extensionThreading),
							fragsDescription.get(c)[3],  rmsCutOff , extensionThreading-overlapFree , extensionThreading+(FRAGSIZE-1) ,-1,overlapFree,3 /*No PG check is needed here*/,
							similarityMatchCO, libs[getPrevLibIndex(c)]);
				libStarts[c] = fragsDescription.get(c)[0] - resStart;
				libEnds[c] = libStarts[c] + (FRAGSIZE-1);
				libManners[c] = -1;
				if (true) { // To print
					System.out.println("FRAGSIZE " + FRAGSIZE);
					for (int d=0; d<Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading, fragsDescription.get(c)[0]+(FRAGSIZE-1)+extensionThreading).length ; d++)
						System.out.print(Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading, fragsDescription.get(c)[0]+(FRAGSIZE-1)+extensionThreading)[d]+ " ");
					System.out.println();
					System.out.println(libStarts[c] + "   SSS   " + libEnds[c]);
					System.out.println("variable Lib Size is:" + libs[c].libSize());
				}
			}
			if (fragsDescription.get(c)[2]==2) {
				int FRAGSIZE = fragsDescription.get(c)[1]-fragsDescription.get(c)[0]+1;
				libs[c] = new LocalFragmentLib(corpus , Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading, fragsDescription.get(c)[1]+extensionThreading), 
							fragsDescription.get(c)[3],  rmsCutOff, extensionThreading , extensionThreading+(FRAGSIZE-1)+overlapFree ,1,overlapFree,3 /*No PG check is needed here*/,
							similarityMatchCO, libs[getPrevLibIndex(c)]);
				libStarts[c] = fragsDescription.get(c)[0] - resStart;
				libEnds[c] = libStarts[c] + (FRAGSIZE-1);
				libManners[c] = 1;
				if (true) { // To print
					System.out.println("FRAGSIZE " + FRAGSIZE);
					for (int d=0; d<Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading, fragsDescription.get(c)[1]+extensionThreading).length ; d++)
						System.out.print(Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading, fragsDescription.get(c)[1]+extensionThreading)[d]+ " ");
					System.out.println();
					System.out.println(libStarts[c] + "   SSS   " + libEnds[c]);
					System.out.println("variable Lib Size is:" + libs[c].libSize());
				}
			}

			
/*    The old version without the extension at the end
			if (fragsDescription.get(c)[2]==-1) {
				int FRAGSIZE = fragsDescription.get(c)[1]-fragsDescription.get(c)[0]+1;
				libs[c] = new LocalFragmentLib(corpus , Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading,fragsDescription.get(c)[0]+(FRAGSIZE-1)) , 
						fragsDescription.get(c)[3],  rmsCutOff, (extensionThreading-overlap) , (extensionThreading+(FRAGSIZE-1)) ,  
						prot,(fragsDescription.get(c)[0]-overlap),overlap,-1,similarityMatchCO);
				libStarts[c] = fragsDescription.get(c)[0] - resStart;
				libEnds[c] = libStarts[c] + (FRAGSIZE-1);
				libManners[c] = -1;
				if (true) { // To print
					System.out.println("FRAGSIZE " + FRAGSIZE);
					for (int d=0; d<Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading,fragsDescription.get(c)[0]+(FRAGSIZE-1)).length ; d++)
						System.out.print(Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading,fragsDescription.get(c)[0]+(FRAGSIZE-1))[d]+ " ");
					System.out.println();
					System.out.println(libStarts[c] + "   SSS   " + libEnds[c]);
					System.out.println("variable Lib Size is:" + libs[c].libSize());
				}
			}
			if (fragsDescription.get(c)[2]==1) {
				int FRAGSIZE = fragsDescription.get(c)[1]-fragsDescription.get(c)[0]+1;
				libs[c] = new LocalFragmentLib(corpus , Protein.getSeqOfProt(prot, fragsDescription.get(c)[1]-(FRAGSIZE-1),fragsDescription.get(c)[1]+extensionThreading) , 
						fragsDescription.get(c)[3],  rmsCutOff, 0 , ((FRAGSIZE-1)+overlap) ,  
						prot,(fragsDescription.get(c)[1]-(FRAGSIZE-1)),overlap,1,similarityMatchCO);
				libStarts[c] = fragsDescription.get(c)[0] - resStart;
				libEnds[c] = libStarts[c] + (FRAGSIZE-1);
				libManners[c] = 1;
				if (true) { // To print
					System.out.println("FRAGSIZE " + FRAGSIZE);
					for (int d=0; d<Protein.getSeqOfProt(prot, fragsDescription.get(c)[1]-(FRAGSIZE-1),fragsDescription.get(c)[1]+extensionThreading).length ; d++)
						System.out.print(Protein.getSeqOfProt(prot, fragsDescription.get(c)[1]-(FRAGSIZE-1),fragsDescription.get(c)[1]+extensionThreading)[d]+ " ");
					System.out.println();
					System.out.println(libStarts[c] + "   SSS   " + libEnds[c]);
					System.out.println("variable Lib Size is:" + libs[c].libSize());
				}
			}
			if (fragsDescription.get(c)[2]==-2) {
				int FRAGSIZE = fragsDescription.get(c)[1]-fragsDescription.get(c)[0]+1;
				libs[c] = new LocalFragmentLib(corpus , Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading, fragsDescription.get(c)[0]+(FRAGSIZE-1)),
						fragsDescription.get(c)[3],  rmsCutOff , extensionThreading-overlap , extensionThreading+(FRAGSIZE-1) ,-1,overlap);
				libStarts[c] = fragsDescription.get(c)[0] - resStart;
				libEnds[c] = libStarts[c] + (FRAGSIZE-1);
				libManners[c] = -1;
				if (true) { // To print
					System.out.println("FRAGSIZE " + FRAGSIZE);
					for (int d=0; d<Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading, fragsDescription.get(c)[0]+(FRAGSIZE-1)).length ; d++)
						System.out.print(Protein.getSeqOfProt(prot, fragsDescription.get(c)[0]-extensionThreading, fragsDescription.get(c)[0]+(FRAGSIZE-1))[d]+ " ");
					System.out.println();
					System.out.println(libStarts[c] + "   SSS   " + libEnds[c]);
					System.out.println("variable Lib Size is:" + libs[c].libSize());
				}
			}
			if (fragsDescription.get(c)[2]==2) {
				int FRAGSIZE = fragsDescription.get(c)[1]-fragsDescription.get(c)[0]+1;
				libs[c] = new LocalFragmentLib(corpus , Protein.getSeqOfProt(prot, fragsDescription.get(c)[0], fragsDescription.get(c)[1]+extensionThreading), 
						fragsDescription.get(c)[3],  rmsCutOff, 0 , (FRAGSIZE-1)+overlap ,1,overlap);
				libStarts[c] = fragsDescription.get(c)[0] - resStart;
				libEnds[c] = libStarts[c] + (FRAGSIZE-1);
				libManners[c] = 1;
				if (true) { // To print
					System.out.println("FRAGSIZE " + FRAGSIZE);
					for (int d=0; d<Protein.getSeqOfProt(prot, fragsDescription.get(c)[0], fragsDescription.get(c)[1]+extensionThreading).length ; d++)
						System.out.print(Protein.getSeqOfProt(prot, fragsDescription.get(c)[0], fragsDescription.get(c)[1]+extensionThreading)[d]+ " ");
					System.out.println();
					System.out.println(libStarts[c] + "   SSS   " + libEnds[c]);
					System.out.println("variable Lib Size is:" + libs[c].libSize());
				}
			}
*/		
		
		}
		// Finally, setting the data straight for the the 'temporaryDistEnergy'
		int lastStart = resStart-1;
		int lastEnd = resEnd + 1;
		for (int c=0 ; c<fragsDescription.size() ; c++) {
			if (libManners[c]==1) {
				startResForClosingPotential[c] = lastStart;
				endResForClosingPotential[c] = resStart + libStarts[c];
				lastEnd = endResForClosingPotential[c];
			}
			if (libManners[c]==-1) {
				endResForClosingPotential[c] = lastEnd;
				startResForClosingPotential[c] = resStart + libEnds[c];
				lastStart = startResForClosingPotential[c];
			}
		}
		
//		DistanceMatrix dm = new DistanceMatrix(prot.atoms(),  2.0, 1.0, 4);
//		pp = RotamericTools.putIntoRot1(prot , dm , rotLib);	
//		try{
//			BufferedWriter bw = new BufferedWriter(new FileWriter("j:/loop_2nd/libData.txt"));
//			for (int c=0 ; c<libs.length ; c++) {
//				if (libManners[c]==-1) {
//					printLibStat(resStart + (libStarts[c]-2),pp,libs[c],0,c,1,bw);
//					printLibStat(resStart + (libStarts[c]-1),pp,libs[c],1,c,2,bw);
//					printLibStat(resStart + (libStarts[c]),pp,libs[c],2,c,3,bw);
//					printLibStat(resStart + (libStarts[c]+1),pp,libs[c],3,c,4,bw);
//				}
//				if (libManners[c]==1) {
//					printLibStat(resStart + (libStarts[c]+3),pp,libs[c],3,c,1,bw);
//					printLibStat(resStart + (libStarts[c]+2),pp,libs[c],2,c,2,bw);
//					printLibStat(resStart + (libStarts[c]+1),pp,libs[c],1,c,3,bw);
//					printLibStat(resStart + (libStarts[c]),pp,libs[c],0,c,4,bw);
//				}
//			}
//			bw.close();
//		}
//		catch(Exception e) {
//			throw new RuntimeException(e.getMessage());
//		}		
//		System.exit(1);
	}

//	private void printLibStat(int resNum, double[][] ppp, LocalFragmentLib lib, int numInFrag,
//			int libNum, int serialNum, BufferedWriter bw) throws Exception {
//		System.out.println(libNum + " " + serialNum + " -99 " + ppp[resNum][0] + " " + ppp[resNum][1] + " 0.0 " + resNum);
//		bw.write(libNum + " " + serialNum + " -99 " + ppp[resNum][0] + " " + ppp[resNum][1] + " 0.0 " + resNum + "\n");		
//		for (int c=0 ; c<200 ; c++) {
//			System.out.println(libNum + " " + serialNum + " " + c + " " + corpus.torsions[lib.libOrig(c)+numInFrag][1] + " " +
//					corpus.torsions[lib.libOrig(c)+numInFrag][2] + " " + lib.libEnergy(c) + " " + resNum);
////			System.out.print(lib.printLib(c));
//			bw.write(libNum + " " + serialNum + " " + c + " " + corpus.torsions[lib.libOrig(c)+numInFrag][1] + " " +
//					corpus.torsions[lib.libOrig(c)+numInFrag][2] + " " + lib.libEnergy(c) + " " + resNum + "\n");
//		}		
//	}
	
	
	
//We need to override this method because it asks for a residue that is +1 behind endRes
protected void updateSeqAndPrePRO() {
	seq = Protein.getSeqOfProt(prot,resStart,resEnd);

	// Finding pre-Prolines in the sequence
	prepro = new boolean[gapLength];
	int[] tmpSEQ = Protein.getSeqOfProt(prot,resStart,resEnd);
	for (int c=0 ; c<gapLength-1 ; c++)
		if ((tmpSEQ[c+1]==PRO) && (tmpSEQ[c]!=GLY) && (tmpSEQ[c]!=PRO))
			prepro[c] = true;
		else 
			prepro[c] = false;
	prepro[gapLength-1] = false;
}

protected int giveFragLength() {
    return -999;
}


@Override
protected boolean continueForNextCall(int libCounter) {
	boolean goodEV = ((evaluateEVextLoop()+evaluateEVinterLoop())<EV_CUTOFF);
	boolean goodTemplateTether = true;
	if (tetherToTemplate>0) {
		for (int c=libStarts[libCounter] ; c<=libEnds[libCounter] ; c++) {
			if (((modelCAs[c].x()-templateCAs[c].x())*(modelCAs[c].x()-templateCAs[c].x()) + 
					 (modelCAs[c].y()-templateCAs[c].y())*(modelCAs[c].y()-templateCAs[c].y()) + 
					 (modelCAs[c].z()-templateCAs[c].z())*(modelCAs[c].z()-templateCAs[c].z())) > tetherToTemplateSquared) {
				goodTemplateTether = false;
				break;
			}			
		}		
	}
	if (goodTemplateTether)
		TETHERaccepted++;
	else
		TETHERrejected++;	
	if (goodEV)
		EVaccepted++;
	else
		EVrejected++;			
	return goodEV && goodTemplateTether;
}

/**
 * This method set a structure of template and model tethered atoms for quick reference. 
 */
private void setTemplateTetherStructure() {
	templateCAs = new Atom[resEnd-resStart+1];
	modelCAs = new Atom[resEnd-resStart+1];
	for (int c=0 ; c<(resEnd-resStart+1) ; c++) {
		templateCAs[c] = protCopy.findAtomInList("CA", resStart+c);
		modelCAs[c] = prot.atoms().findAtomInList("CA", resStart+c);
	}
}

public void setTetherToTemplate(double newTether) {
	System.out.println("Tethereing to template with a cutoff of " + newTether + " Angs.");
	tetherToTemplate = newTether;
	tetherToTemplateSquared = newTether*newTether;
}

private int getPrevLibIndex(int c) { // 4/26/2010 - this is a patch that will allow for the testing of the prevLib concept.
	return c-2;
}

}

