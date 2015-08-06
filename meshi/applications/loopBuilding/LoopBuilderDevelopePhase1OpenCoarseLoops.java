package meshi.applications.loopBuilding;

import java.util.Vector;

import meshi.applications.corpus.Corpus;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.ROT1solvation.CBSolvationCreator;
import meshi.energy.ROT1solvation.ROT1SolvationCreator;
import meshi.energy.ROT1solvation.ROT1SolvationEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;

public class LoopBuilderDevelopePhase1OpenCoarseLoops extends
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

	public LoopBuilderDevelopePhase1OpenCoarseLoops(CommandList commands,
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
			allResults.add(new BasicLoopResult(results.evEnergy,results.propEnergy,results.bbHBenergy,0.0/*score is ignored*/,results.RMS,fragRank,saveLoopCoordinates(),saveLoopCoordinatesNCaC(),savePhiPsiOfLoop()));
		}
	}
	
	public void developementAnalysis(int moduleBase) {
		EnergyCreator[] energyCreators = {new ROT1SolvationCreator(1.0,"4.0",false),
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
		energy = new TotalEnergy(prot, new DistanceMatrix(prot.atoms(),  9.0, 0.1, 4), energyCreators, commands);
		rot1_040 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[0]);
		rot1_045 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[1]);
		rot1_050 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[2]);
		rot1_055 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[3]);
		rot1_060 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[4]);
		rot1_065 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[5]);
		rot1_070 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[6]);
		rot1_080 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[7]);
		rot1_090 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[8]);
		rot1_100 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[9]);
		rot1_120 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[10]);
		rot1_140 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[11]);
		cb_040 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[12]);
		cb_045 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[13]);
		cb_050 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[14]);
		cb_055 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[15]);
		cb_060 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[16]);
		cb_065 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[17]);
		cb_070 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[18]);
		cb_080 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[19]);
		cb_090 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[20]);
		cb_100 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[21]);
		cb_120 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[22]);
		cb_140 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[23]);	

		for (int modelNum=0 ; modelNum<allResults.size() ; modelNum++) {
			restoreLoopCoordinates(allResults.get(modelNum).coors);
			try {
				energy.update();
			}
			catch (Exception e) {
				System.out.println("Warning: an error occured in the energy calculation of the development stage. The loop will be dicarded.\n the error is:" + e.getMessage());
			}
			// The data on the loop
			System.out.print("777777 " + (moduleBase+modelNum) + " " +
					fmt2.format(allResults.get(modelNum).rms) + " " +
					fmt2.format(allResults.get(modelNum).evEnergy) + " " +
					fmt2.format(allResults.get(modelNum).propEnergy) + " " + 
					fmt2.format(allResults.get(modelNum).bbHBEnergy) + "          " +
					fmt2.format(rot1_040.evaluate()) + " " + 
					fmt2.format(rot1_045.evaluate()) + " " + 
					fmt2.format(rot1_050.evaluate()) + " " + 
					fmt2.format(rot1_055.evaluate()) + " " + 
					fmt2.format(rot1_060.evaluate()) + " " + 
					fmt2.format(rot1_065.evaluate()) + " " + 
					fmt2.format(rot1_070.evaluate()) + " " + 
					fmt2.format(rot1_080.evaluate()) + " " + 
					fmt2.format(rot1_090.evaluate()) + " " + 
					fmt2.format(rot1_100.evaluate()) + " " + 
					fmt2.format(rot1_120.evaluate()) + " " + 
					fmt2.format(rot1_140.evaluate()) + "         " + 
					fmt2.format(cb_040.evaluate()) + " " + 
					fmt2.format(cb_045.evaluate()) + " " + 
					fmt2.format(cb_050.evaluate()) + " " + 
					fmt2.format(cb_055.evaluate()) + " " + 
					fmt2.format(cb_060.evaluate()) + " " + 
					fmt2.format(cb_065.evaluate()) + " " + 
					fmt2.format(cb_070.evaluate()) + " " + 
					fmt2.format(cb_080.evaluate()) + " " + 
					fmt2.format(cb_090.evaluate()) + " " + 
					fmt2.format(cb_100.evaluate()) + " " + 
					fmt2.format(cb_120.evaluate()) + " " + 
					fmt2.format(cb_140.evaluate()) + "         ");
			System.out.print(fmt2.format(rot1_040.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(rot1_045.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(rot1_050.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(rot1_055.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(rot1_060.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(rot1_065.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(rot1_070.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(rot1_080.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(rot1_090.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(rot1_100.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(rot1_120.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(rot1_140.getEnergyHydrophobicWeightedSA()) + "         " + 
					fmt2.format(cb_040.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(cb_045.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(cb_050.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(cb_055.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(cb_060.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(cb_065.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(cb_070.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(cb_080.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(cb_090.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(cb_100.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(cb_120.getEnergyHydrophobicWeightedSA()) + " " + 
					fmt2.format(cb_140.getEnergyHydrophobicWeightedSA()) + "         ");
			System.out.print(fmt2.format(rot1_040.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(rot1_045.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(rot1_050.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(rot1_055.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(rot1_060.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(rot1_065.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(rot1_070.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(rot1_080.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(rot1_090.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(rot1_100.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(rot1_120.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(rot1_140.getEnergyPolarWeightedSA()) + "         " + 
					fmt2.format(cb_040.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(cb_045.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(cb_050.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(cb_055.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(cb_060.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(cb_065.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(cb_070.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(cb_080.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(cb_090.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(cb_100.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(cb_120.getEnergyPolarWeightedSA()) + " " + 
					fmt2.format(cb_140.getEnergyPolarWeightedSA()) + "         ");
			System.out.print(fmt2.format(rot1_040.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(rot1_045.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(rot1_050.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(rot1_055.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(rot1_060.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(rot1_065.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(rot1_070.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(rot1_080.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(rot1_090.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(rot1_100.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(rot1_120.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(rot1_140.getEnergyRepresentativePolar()) + "         " + 
					fmt2.format(cb_040.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(cb_045.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(cb_050.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(cb_055.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(cb_060.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(cb_065.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(cb_070.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(cb_080.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(cb_090.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(cb_100.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(cb_120.getEnergyRepresentativePolar()) + " " + 
					fmt2.format(cb_140.getEnergyRepresentativePolar()) + "         ");
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
	}
	
	public void setProt(Protein newProt) {
		prot = newProt;
	}

}
