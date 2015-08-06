package meshi.applications.loopBuilding;

import java.util.Calendar;
import java.util.Vector;

import meshi.applications.corpus.Corpus;
import meshi.molecularElements.Atom;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;

public class LoopBuilderUserDefinedFragmentsStochasticWithClosure extends LoopBuilderUserDefinedFragmentsStochastic {
	
	protected final double GAP_CUTOFF=1.0; // The building tree is pruned if the gap score is above this value 	
	
	private double closureToleranceSquared = -99999.0;  // We precalculate the square to save computaion time
	protected Atom atom1C, atom2N; // The closing atoms over the final peptide bond
	
	// Stats
	private int numberReachClosure = 0;
	private int DISaccepted = 0;
	private int DISrejected = 0;

	public LoopBuilderUserDefinedFragmentsStochasticWithClosure(CommandList commands,
			String writePath, Corpus corpus, Protein prot, Protein ref,
			int resStart, int resEnd, double rmsMatchCO, double rmsCutOff,
			Vector<int[]> fragsDescription, double closureTolerance) {
		super(commands, writePath, corpus, prot, ref, resStart, resEnd,
				rmsMatchCO, rmsCutOff, fragsDescription);
		this.closureToleranceSquared = closureTolerance*closureTolerance;
		System.out.println("Closure tolerance: " + closureTolerance);
	}

	protected void setClosureAtoms() {
		int lastInd = libs.length - 1;
		int closureRes;
		if (libManners[lastInd]==-1) // Ends from the N direction
			closureRes = resStart + libEnds[lastInd];
		else   // Ends from the C direction
			closureRes = resStart + libStarts[lastInd]-1;
		atom1C = prot.atoms().findAtomInList("C", closureRes);
		atom2N = prot.atoms().findAtomInList("N", closureRes+1);
		System.out.println("Closure atoms:\n"+atom1C+"\n"+atom2N);
		System.out.print("Closure residues:\n"+atom1C.residueNumber()+" 9999999999\n"+atom2N.residueNumber()+" 9999999999\n");		
	}

	
	private final double[][] tolerances = {
			{2.5 , 60.0*Math.PI/180.0}, //0
			{2.75 , 60.0*Math.PI/180.0}, //1 
			{3.0 , 60.0*Math.PI/180.0}, //2
			{3.25 , 60.0*Math.PI/180.0}, //3
			{3.5 , 70.0*Math.PI/180.0}, //4
			{3.6 , 90.0*Math.PI/180.0}, //5
			{3.7 , 120.0*Math.PI/180.0}, //6
			{3.8 , 999.0*Math.PI/180.0}, //7
			{3.9 , 999.0*Math.PI/180.0}, //8
			{4.0 , 999.0*Math.PI/180.0}, //9
			{4.0 , 999.0*Math.PI/180.0}, //10
			{4.0 , 999.0*Math.PI/180.0}, //11
	};

	protected void buildCoarsFragments() {
		long lastCheckPointTime = Calendar.getInstance().getTimeInMillis();
		int numberOfModelsLastCheckPoint = 0;
		long startingTime = Calendar.getInstance().getTimeInMillis();
		long nowTime = startingTime;
		int toleranceCounter = 0;
		while (!gotEnoughCoarseModels && ((nowTime-startingTime)/1000 < MAX_SECONDS_TO_RUN)) {
			addFromLib(0);
			nowTime = Calendar.getInstance().getTimeInMillis();
			if ((nowTime-lastCheckPointTime)/1000 > (MAX_SECONDS_TO_RUN/10.0)) { // Checkpoint
				if ((allResults.size()-numberOfModelsLastCheckPoint) < (maxNumberOfLoopGenerated/10)) { // Progress NOT satisfactory
					toleranceCounter++;
					closureToleranceSquared = tolerances[toleranceCounter][0]*tolerances[toleranceCounter][0];
					similarityMatchCO = tolerances[toleranceCounter][1];
					System.out.println("\n\n\nThe closure tolerance was upped to : " + Math.sqrt(closureToleranceSquared) + "\n" +
							"The stalk similarity tolerance was upped to : " + (similarityMatchCO/Math.PI*180.0) + "\n\n\n");				
				}
				lastCheckPointTime = nowTime;
				numberOfModelsLastCheckPoint = allResults.size();				
			}
		}
		printStats();
	}	
	
	/* I used this code for a long time until 13/3/2010, but then changed it for the more refined version above.
	protected void buildCoarsFragments() {
		long startingTime = Calendar.getInstance().getTimeInMillis();
		long nowTime = startingTime;
		while (!gotEnoughCoarseModels && ((nowTime-startingTime)/1000 < MAX_SECONDS_TO_RUN)) {
			addFromLib(0);
			nowTime = Calendar.getInstance().getTimeInMillis();
		}
		if (allResults.size()<(maxNumberOfLoopGenerated/5)) {
			closureToleranceSquared = (Math.sqrt(closureToleranceSquared)+1.0)*(Math.sqrt(closureToleranceSquared)+1.0);
			System.out.println("\n\n\nThe closure tolerance was upped to : " + Math.sqrt(closureToleranceSquared) + "\n\n\n");			
		}
		startingTime = Calendar.getInstance().getTimeInMillis();
		nowTime = startingTime;
		while (!gotEnoughCoarseModels && ((nowTime-startingTime)/1000 < MAX_SECONDS_TO_RUN)) {
			addFromLib(0);
			nowTime = Calendar.getInstance().getTimeInMillis();
		}
		printStats();
	}	
	*/
	
	/* I used this overide to get 5000 loops on the Honig closure set 
	protected void buildCoarsFragments() {
		int upClosureCounter = 0;
		for (; (upClosureCounter<3) && !gotEnoughCoarseModels ; upClosureCounter++) {
			long startingTime = Calendar.getInstance().getTimeInMillis();
			long nowTime = startingTime;
			while (!gotEnoughCoarseModels && ((nowTime-startingTime)/1000 < MAX_SECONDS_TO_RUN)) {
				loopCompleted = false;
				addFromLib(0);
				nowTime = Calendar.getInstance().getTimeInMillis();
			}
			if (!gotEnoughCoarseModels) { 
				closureToleranceSquared = (Math.sqrt(closureToleranceSquared)+0.5)*(Math.sqrt(closureToleranceSquared)+0.5);
				System.out.println("\n\n\nThe closure tolerance was upped to : " + Math.sqrt(closureToleranceSquared) + "\n\n\n");			
			}
		}
		if (upClosureCounter>2) 
			System.out.println("\n\n\nWARNING: the maximal number of tolerance increments reached.\n\n\n");
		System.out.println("\n\n\nThe final closure tolerance: " + Math.sqrt(closureToleranceSquared) + "\n\n\n");
		printStats();
	}
*/
	
	protected boolean loopClosingCondition() {
		numberReachClosure++;
		// First: the N and C must be within the tolerance 
		if (((atom1C.x()-atom2N.x())*(atom1C.x()-atom2N.x()) + 
			 (atom1C.y()-atom2N.y())*(atom1C.y()-atom2N.y()) + 
			 (atom1C.z()-atom2N.z())*(atom1C.z()-atom2N.z())) < closureToleranceSquared) {
					return true;
		}
		return false;
	}

	@Override
	protected boolean continueForNextCall(int libCounter) {
		boolean goodDis = (temporaryDistEnergy(startResForClosingPotential[libCounter],endResForClosingPotential[libCounter])<GAP_CUTOFF);
		if (goodDis)
			DISaccepted++;
		else
			DISrejected++;
		return (goodDis && super.continueForNextCall(libCounter));
	}

	protected void 	printStats() {
	    System.out.println("Number of closure assessed: " + numberReachClosure);	
	    System.out.println("Number of distance acceptances: " + DISaccepted);	
	    System.out.println("Number of distance rejectences: " + DISrejected);	
	    System.out.println("Number of EV acceptances: " + EVaccepted);	
	    System.out.println("Number of EV rejectences: " + EVrejected);	
	    System.out.println("Number of TETHER acceptances: " + TETHERaccepted);	
	    System.out.println("Number of TETHER rejectences: " + TETHERrejected);		    
	}
		
}
