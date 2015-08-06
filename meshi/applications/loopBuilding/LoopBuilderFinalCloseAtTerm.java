package meshi.applications.loopBuilding;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.Vector;

import meshi.applications.corpus.Corpus;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.ROT1solvation.CBSolvationCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;

public class LoopBuilderFinalCloseAtTerm extends
LoopBuilderUserDefinedFragmentsStochasticWithClosure {
	
	protected int loopType = -999;
	
	public LoopBuilderFinalCloseAtTerm(CommandList commands,
			String writePath, Corpus corpus, Protein prot, Protein ref,
			int resStart, int resEnd, double rmsMatchCO, double rmsCutOff,
			Vector<int[]> fragsDescription, double closureTolerance, int loopType) {
		super(commands, writePath, corpus, prot, ref, resStart, resEnd,
				rmsMatchCO, rmsCutOff, fragsDescription, closureTolerance);
		this.loopType = loopType;
	}

	/** 
	 * If the loop is constrained in both ends we run the 'super' method in 'LoopBuilderUserDefinedFragmentsStochasticWithClosure'.
	 * However, if this is a terminal loop we set the closure atoms to be at the loop stalk. This guarantee they are always accepted.
	 */
	protected void setClosureAtoms() {
		if (loopType==0)
			super.setClosureAtoms();
		else {
			if (loopType==-1) { // The loop is constrained in the Nterm
				atom1C = prot.atoms().findAtomInList("C", resEnd);
				atom2N = prot.atoms().findAtomInList("N", resEnd+1);
			}
			else { // The loop is constrained in the Cterm
				atom1C = prot.atoms().findAtomInList("C", resStart-1);
				atom2N = prot.atoms().findAtomInList("N", resStart);				
			}
			System.out.println("Closure atoms:\n"+atom1C+"\n"+atom2N);
			System.out.print("Closure residues:\n"+atom1C.residueNumber()+" 9999999999\n"+atom2N.residueNumber()+" 9999999999\n");
		}
	}
	
	// Build the local frag libraries.
	protected void buildLibraries() {
		super.buildLibraries();
		// Setting the data straight for the the 'temporaryDistEnergy'
		if (loopType==-1) { // The loop is constrained in the Nterm
			for (int c=0 ; c<fragsDescription.size() ; c++) {
				startResForClosingPotential[c] = resStart+libEnds[c];
				endResForClosingPotential[c] = resEnd+1;
			}			
		}
		else if (loopType==1) { // The loop is constrained in the Cterm
			for (int c=0 ; c<fragsDescription.size() ; c++) {
				startResForClosingPotential[c] = resStart-1;
				endResForClosingPotential[c] = resStart+libStarts[c];
			}						
		}
	}	
	
	
	protected void analyzedCloseLoop() {
		InternalLoopData results;
		try {
			results = evaluateScores(resStart ,resEnd); 	
		}
		catch (Exception e) {
			System.out.println("Warning: an error occured in the energy calculation. The loop will be dicarded.\n the error is:" + e.getMessage());
			return;
		}
		if (results.evEnergy<EV_CUTOFF) {
			numberOfCoarseGenerated++;
			if (numberOfCoarseGenerated>maxNumberOfLoopGenerated)
				gotEnoughCoarseModels = true;
			if (numberOfCoarseGenerated/howOftenToPrintFinalLoop == ((int) numberOfCoarseGenerated/(1.0*howOftenToPrintFinalLoop))) {
				System.out.print(PRINTING_HEADER_COARSE + " " + numberOfCoarseGenerated + " " +
						fmt2.format(results.RMS) + " " +
						fmt2.format(results.propEnergy) + " " +
						fmt2.format(results.evEnergy) + " " + 
						fmt2.format(results.bbHBenergy) + "         ");
				for (int tmpc=0 ; tmpc<fragRank.length ; tmpc++) {
					System.out.print(fragRank[tmpc]+" ");
					//score += 0.0*(Math.log(1+fragRank[tmpc]));
				}
				System.out.print("       ");
				for (int tmpc=0 ; tmpc<actuallyTaken.length ; tmpc++) {
					System.out.print(actuallyTaken[tmpc]+" ");
				}
				System.out.println();
			}
			allResults.add(new BasicLoopResult(results.evEnergy,results.propEnergy,results.bbHBenergy,0.0/*score is ignored*/,results.RMS,fragRank,saveLoopCoordinates(),saveLoopCoordinatesNCaC(),savePhiPsiOfLoop()));
		}
	}	
	
	public void cbAnalysis(int numToPrint) {
		EnergyCreator[] energyCreators = {new CBSolvationCreator(1.0,"9.0",false)};
		energy = new TotalEnergy(prot, new DistanceMatrix(prot.atoms(),  9.0, 0.1, 4), energyCreators, commands);
		double cbEne  = 0.0;

		double[] scores = new double[allResults.size()]; // For the ranking of the total scores
		for (int modelNum=0 ; modelNum<allResults.size() ; modelNum++) {
			restoreLoopCoordinates(allResults.get(modelNum).coors);
			cbEne = energy.evaluate();
			allResults.get(modelNum).setScore(cbEne);
			scores[modelNum] = 0.4*allResults.get(modelNum).evEnergy +
			0.4*allResults.get(modelNum).propEnergy +
			1.0*allResults.get(modelNum).bbHBEnergy +
			1.5*cbEne;
		}
		int numberToSave = Math.min(allResults.size(), numToPrint);
		int[] sortedScores = findTopMinArray(scores, numberToSave, infiniteEnergy);
		
		// Printing to disk the best ones		
		System.out.println("Printing the final " + numberToSave + "models:\n----------------------------------------\n");
		DecimalFormat fmt = new DecimalFormat("0.##");    
		for (int modelNum=0 ; modelNum<numberToSave ; modelNum++) {
			System.out.println(modelNum + " " + sortedScores[modelNum] + " " + 
					fmt.format(scores[sortedScores[modelNum]]) + " " + 
					fmt.format(allResults.get(sortedScores[modelNum]).score) + " " + 
					fmt.format(allResults.get(sortedScores[modelNum]).bbHBEnergy) + " " + 
					fmt.format(allResults.get(sortedScores[modelNum]).propEnergy) + " " + 
					fmt.format(allResults.get(sortedScores[modelNum]).evEnergy) + " " + 
					fmt.format(allResults.get(sortedScores[modelNum]).rms));
			restoreLoopCoordinates(allResults.get(sortedScores[modelNum]).coors);
			// Printing the coordinates
			try{
				BufferedWriter bw = new BufferedWriter(new FileWriter(writePath+"/phase1/"+modelNum+".pdb"));
				for (int res=resStart ; res<=resEnd ; res++) {
					for (int atInd=0 ; atInd<prot.residue(res).atoms().size() ; atInd++) {
						bw.write(prot.residue(res).atoms().atomAt(atInd).toString() + "\n");
					}
				}
				bw.close();
			}
			catch(Exception e) {
				throw new RuntimeException(e.getMessage());
			}
			// Printing the data
			try{
				BufferedWriter bw = new BufferedWriter(new FileWriter(writePath+"/phase1/"+modelNum+".data"));
				for (int res=resStart ; res<=resEnd ; res++) {
					bw.write(res + " " + 
							fmt.format(allResults.get(modelNum).myPhiPsi[res-resStart][0]) + " " + 
							fmt.format(allResults.get(modelNum).myPhiPsi[res-resStart][1]) + "\n");
				}
				bw.close();
			}
			catch(Exception e) {
				throw new RuntimeException(e.getMessage());
			}
		}
	}
	
	public void setProt(Protein newProt) {
		prot = newProt;
	}

}
