package meshi.applications.loopBuilding;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.Vector;

import meshi.applications.corpus.Corpus;
import meshi.energy.ROT1solvation.ROT1SolvationEnergy;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;

public class LoopBuilderDevelopePhase2OpenCoarseLoops extends
LoopBuilderUserDefinedFragmentsStochasticWithClosure {
	
	protected ROT1SolvationEnergy rot1_040 = null;
	protected ROT1SolvationEnergy rot1_045 = null;
	protected ROT1SolvationEnergy rot1_050 = null;
	protected ROT1SolvationEnergy rot1_055 = null;
	protected ROT1SolvationEnergy rot1_060 = null;
	protected ROT1SolvationEnergy rot1_065 = null;
	protected ROT1SolvationEnergy rot1_070 = null;
	protected ROT1SolvationEnergy rot1_080 = null;
	protected ROT1SolvationEnergy rot1_090 = null;
	protected ROT1SolvationEnergy rot1_100 = null;
	protected ROT1SolvationEnergy rot1_120 = null;
	protected ROT1SolvationEnergy rot1_140 = null;
	protected ROT1SolvationEnergy cb_040 = null;
	protected ROT1SolvationEnergy cb_045 = null;
	protected ROT1SolvationEnergy cb_050 = null;
	protected ROT1SolvationEnergy cb_055 = null;
	protected ROT1SolvationEnergy cb_060 = null;
	protected ROT1SolvationEnergy cb_065 = null;
	protected ROT1SolvationEnergy cb_070 = null;
	protected ROT1SolvationEnergy cb_080 = null;
	protected ROT1SolvationEnergy cb_090 = null;
	protected ROT1SolvationEnergy cb_100 = null;
	protected ROT1SolvationEnergy cb_120 = null;
	protected ROT1SolvationEnergy cb_140 = null;

	public LoopBuilderDevelopePhase2OpenCoarseLoops(CommandList commands,
			String writePath, Corpus corpus, Protein prot, Protein ref,
			int resStart, int resEnd, double rmsMatchCO, double rmsCutOff,
			Vector<int[]> fragsDescription, double closureTolerance) {
		super(commands, writePath, corpus, prot, ref, resStart, resEnd,
				rmsMatchCO, rmsCutOff, fragsDescription, closureTolerance);
	}
/*
	protected EnergyCreator[] getEnergyCreatorArray() {
		EnergyCreator[] energyCreators = {new  SoftExcludedVolCreator(1.0 , 9 , 1.0),
				new  SoftExcludedVolCreator(1.0 , 10 , 1.0),
				new SimpleHydrogenBond_Dahiyat_LowAccuracy_BBonly_Creator(1.0),
				new ROT1SolvationCreator(1.0,"4.0",false),
				new ROT1SolvationCreator(1.0,"4.5",false),
				new ROT1SolvationCreator(1.0,"5.0",false),
				new ROT1SolvationCreator(1.0,"5.5",false),
				new ROT1SolvationCreator(1.0,"6.0",false),
				new ROT1SolvationCreator(1.0,"6.5",false),
				new ROT1SolvationCreator(1.0,"7.0",false),
				new ROT1SolvationCreator(1.0,"8.0",false),
				new ROT1SolvationCreator(1.0,"9.0",false),
				new ROT1SolvationCreator(1.0,"10.0",false),
				new ROT1SolvationCreator(1.0,"12.0",false),
				new ROT1SolvationCreator(1.0,"14.0",false),
				new CBSolvationCreator(1.0,"4.0",false),
				new CBSolvationCreator(1.0,"4.5",false),
				new CBSolvationCreator(1.0,"5.0",false),
				new CBSolvationCreator(1.0,"5.5",false),
				new CBSolvationCreator(1.0,"6.0",false),
				new CBSolvationCreator(1.0,"6.5",false),
				new CBSolvationCreator(1.0,"7.0",false),
				new CBSolvationCreator(1.0,"8.0",false),
				new CBSolvationCreator(1.0,"9.0",false),
				new CBSolvationCreator(1.0,"10.0",false),
				new CBSolvationCreator(1.0,"12.0",false),
				new CBSolvationCreator(1.0,"14.0",false)};
		return energyCreators;
	}	

	protected void resetTotalEnergy() {
		super.resetTotalEnergy();
	}
*/
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
			if (numberOfCoarseGenerated/10 == numberOfCoarseGenerated/10.0) {
				allResults.add(new BasicLoopResult(results.evEnergy,results.propEnergy,results.bbHBenergy,0.0/*score is ignored*/,results.RMS,fragRank,saveLoopCoordinates(),saveLoopCoordinatesNCaC(),savePhiPsiOfLoop()));
			}
		}
	}
	
	public void developementAnalysis(int moduleBase) {
		DecimalFormat fmt = new DecimalFormat("0.##");    
		for (int modelNum=0 ; modelNum<allResults.size() ; modelNum++) {
			restoreLoopCoordinates(allResults.get(modelNum).coors);
			// Printing the coordinates
			try{
				BufferedWriter bw = new BufferedWriter(new FileWriter(writePath+"/"+(moduleBase+modelNum)+".pdb"));
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
				BufferedWriter bw = new BufferedWriter(new FileWriter(writePath+"/"+(moduleBase+modelNum)+".data"));
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
	
	public BasicLoopResultVector getAllResults() {
		return allResults;
	}

}
