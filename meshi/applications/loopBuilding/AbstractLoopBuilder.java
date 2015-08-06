package meshi.applications.loopBuilding;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Vector;

import meshi.applications.corpus.Corpus;
import meshi.applications.localFragLibs.LocalFragmentLib;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.compositePropensity2DwithPP.CompositePropensity2DEnergy;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBondEnergy;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_LowAccuracy_BBonly_Creator;
import meshi.energy.softExcludedVol.SoftExcludedVol;
import meshi.energy.softExcludedVol.SoftExcludedVolCreator;
import meshi.energy.softExcludedVol.SoftExcludedVolEnergyElement;
import meshi.energy.softExcludedVol.SoftExcludedVolParameters;
import meshi.energy.softExcludedVol.SoftExcludedVolParametersList;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.parameters.AtomTypes;
import meshi.parameters.MeshiPotential;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.UpdateableException;
import meshi.util.file.MeshiWriter;
import meshi.util.rotamericTools.RotamericTools;
import programs.SCMOD;

/**
 *This class build a loop. More details later...
 **/

public abstract class AbstractLoopBuilder implements Residues,AtomTypes,CompositeTorsionsDefinitions,MeshiPotential,
KeyWords { 

//	User defined hard-coded parameters:
	protected final int extensionThreading = 9 /*3*/; // must be LARGER than overlap
	protected final int overlapStalk = 2; 
	protected final int overlapFree = 1; 
	protected final double DEFROST_DIST = 8.0; // All residues whose Cb is closer than this to the loop, will be subjected to side-chain modeling.	

	protected final double infiniteEnergy = 1e9; 	
	protected int maxNumberOfLoopGenerated = -999; // The maximal number of coarse loop to be generated. Once this number is reached, farther conformational search ends.
	protected int maxNumberOfLoopToPrint = -999;  // Number of Coarse loops to write to disk.
	protected int takeForSCMOD = -999; // How many coarse models will be chosen for further SCMOD analysis
	protected Corpus corpus = null;
	protected double similarityMatchCO = -999;
	protected int gapLength = -1;  // As the name implies
	protected int resStart = -1; // [resStart,resEnd] are MISSING in the protein
	protected int resEnd = -1; // [resStart,resEnd] are MISSING in the protein
	protected Protein prot=null; // This is the instance that is modeled
	protected Protein ref=null;  // If this is available, than RMS can be calculated
	protected AtomList protCopy = null; // A copy of the original positions of the 
	protected DunbrackLib rotLib = null;
	protected String writePath = "";
	protected double[] energies = null;
	protected double[][] refCoors = null;  // The coordinates of the N,Ca,C loop atoms in the native reference (if available)
	protected double[][] pp = null;   // The {Phi,Psi} of all the residues in the protein
	protected DecimalFormat fmt2 = new DecimalFormat("0.###");
	protected int[] seq;     // Sequence of the missing residues in the loop
	protected boolean prepro[] = null;  // True means the residue is pre-proline


//	Energy related fields
	protected CommandList commands;
	protected TotalEnergy energy = null;
	protected SoftExcludedVol energyTermEVloopFinal = null;
	protected SoftExcludedVol energyTermEVloopAssembly = null;
	protected CompositePropensity2DEnergy energyTermProp = null;
	protected SimpleHydrogenBondEnergy energyTermHB = null;
	private Vector<LocalEVdata> evData_interLoop;
	private Vector<LocalEVdata> evData_extLoop;
	private Vector<LocalEVdata> touchData_extLoop;


//	The object that stores all the loops built so far, and the total number of coarse loops generated.
	protected BasicLoopResultVector allResults = null;
	int numberOfCoarseGenerated = 0;
	protected boolean loopCompleted = false;
	protected boolean gotEnoughCoarseModels = false;

//	Libraries, and their limits. See the first lines in the 'buildLibraries' methods
	protected LocalFragmentLib[] libs = null;  
	protected int[] libStarts = null;  // The indices in the loop (first missing res is 0) where the frags (not overlaps) of each lib strat.
	protected int[] libEnds = null;    // The same as above, but where they end
	protected int[] libManners = null;  // The manner of each lib. (-1) - overlapping in the N terminus. (1) - overlapping in the C terminus.
	protected int[] fragRank = null;   // When a loop is closed, 'fragRank' gives the indices in the libraries of the composing loop fragments.
	protected int[] actuallyTaken = null;  // How deep in the library has the search gone. 
	protected int[] startResForClosingPotential = null; // For each lib, this and the next array can be used as the parameters for the term 'temporaryDistEnergy'
	protected int[] endResForClosingPotential = null;
	
//	Auxiliary variables for the structural alignment
	protected double[][] forAlign1=null;
	protected double[][] forAlign2=null;

	/**
	 * @param commands - obvious
	 * @param writePath - The output files will be written to here as: 0.pdb , 1.pdb, 2.pdb ...
	 * @param corpus - A corpus instance to build frag libs from.
	 * @param prot - The loops will be built on this protein, which must have instances for all the atoms in the loop.
	 * @param ref - This is optional. If it is not 'null' then global-RMS calculations are done in various stages.
	 * @param resStart - [resStart,resEnd] are MISSING in the protein
	 * @param resEnd -   " " " " " " " " " " " " " " " " " " " " " " 
	 * @param similarityMatchCO - When choosing a fragment, its overlap should be better then this value with the existing scaffold.
	 */

	public AbstractLoopBuilder(CommandList commands, String writePath, Corpus corpus, Protein prot, Protein ref, 
			int resStart, int resEnd, double similarityMatchCO) {  
		this.commands = commands;
		this.writePath = writePath;
		this.resStart = resStart;
		this.resEnd = resEnd;	
		this.prot = prot;
		this.ref = ref;
		this.corpus = corpus;
		this.similarityMatchCO = Math.PI/180.0*similarityMatchCO;
		protCopy = prot.atoms().duplicate();
		rotLib = new DunbrackLib(commands, 1.0 , 2);
		gapLength = resEnd-resStart+1; 
	}

	/**
	 * Build an ensamble of coarse loops. 
	 * @param maxNumberOfLoopGenerated - After this many final loop are generated, the method returns.
	 * @param maxNumberOfLoopToPrint - These first loops in 'allResults' will be printed.
	 */
	public void runCoarse(int maxNumberOfLoopGenerated, int maxNumberOfLoopToPrint) {
		this.maxNumberOfLoopGenerated = maxNumberOfLoopGenerated;
		this.maxNumberOfLoopToPrint = maxNumberOfLoopToPrint;
		allResults = new BasicLoopResultVector();
		numberOfCoarseGenerated = 0;
		gotEnoughCoarseModels = false;

		// Finding the sequence
		updateSeqAndPrePRO();

		// Building the local frag libraries
		buildLibraries();
//		System.exit(0);

		// Setting the atoms for the closure
		setClosureAtoms();

		// Building the energy element:
		// ----------------------------
		// Put every residue in Rot1:
		prot.defrost();
		DistanceMatrix dm = new DistanceMatrix(prot.atoms(),  2.0, 1.0, 4);
		pp = RotamericTools.putIntoRot1(prot , dm , rotLib);
		dm = null;
		// Freezing every thing but the loop:
		prot.freeze();
		for (int c=resStart ; c<=resEnd ; c++)
			prot.residue(c).atoms().defrost();
		resetTotalEnergy();


		// For all RMS calculations against the native structure (if available)
		if (ref != null) {
			updateRefCoors();
		}

		// Putting the current loop atoms away so they don't bother the fragments being built
		putInitialLoopFar();    

		buildCoarsFragments();

		// Printing the Coarse models.
		System.out.println("\nPrinting the coarse models...\n");
		int numberToSave = Math.min(allResults.size(), maxNumberOfLoopToPrint);
		for (int modelNum=0 ; modelNum<numberToSave ; modelNum++) {
			// Restoring the loop coordinates
			restoreLoopCoordinates(allResults.get(modelNum).coors);
			// Printing to file
			AtomList loopAtoms = new AtomList();
			for (int res=resStart ; res<=resEnd ; res++)
				loopAtoms.add(prot.residue(res).atoms());
			try {
				loopAtoms.print(new MeshiWriter(writePath+"/coarse/"+modelNum+".pdb"));
			} catch (IOException e) {
				System.out.println("Error writing:" + writePath+"/coarse/"+modelNum+".pdb");
			}
			// Printing to the screen
			System.out.println("555555 " + modelNum + " " + allResults.get(modelNum).score + " " + allResults.get(modelNum).evEnergy + " " + allResults.get(modelNum).calcRMS(allResults.get(0)));
		}
	}

	protected void printFullCoarseModels() {
		// Printing the Coarse models.
		System.out.println("\nPrinting the coarse models...\n");
		int numberToSave = Math.min(allResults.size(), maxNumberOfLoopToPrint);
		for (int modelNum=0 ; modelNum<numberToSave ; modelNum++) {
			// Restoring the loop coordinates
			restoreLoopCoordinates(allResults.get(modelNum).coors);
			// Printing to file
			try {
				prot.atoms().print(new MeshiWriter(writePath+"/coarse/"+modelNum+".pdb"));
			} catch (IOException e) {
				System.out.println("Error writing:" + writePath+"/coarse/"+modelNum+".pdb");
			}
			// Printing to the screen
			System.out.println("999999 " + modelNum + " " + allResults.get(modelNum).score + " " + allResults.get(modelNum).evEnergy + " " + allResults.get(modelNum).calcRMS(allResults.get(0)));
		}
	}

	
	/**
	 * Processing the coarse models and running SCMOD
	 * ----------------------------------------------
	 */
	public void runSCMOD(int takeForSCMOD, double EV_CUTOFF) {
		this.takeForSCMOD = takeForSCMOD;
		if (takeForSCMOD>maxNumberOfLoopToPrint) {
			throw new RuntimeException("You asked for more SCMOD models than are saved (as defined in the toPrint)");
		}

		// First: removing EVs above median   
		double[] EVs = new double[allResults.size()]; // For finding the median of the EV
		for (int c=0; c<EVs.length ; c++)
			EVs[c] = allResults.get(c).evEnergy;
		int[] sortedEVs = findTopMinArray(EVs, EVs.length, infiniteEnergy);
		double medianEV = EVs[sortedEVs[sortedEVs.length/2]];
		System.out.println("The median EV of the loops is: " + medianEV);
		if (medianEV<1.5*EV_CUTOFF)
			medianEV = 1.5*EV_CUTOFF;

		// Second: Sorting the scores
		double largeEnergyAddition = 100000;
		energies = new double[allResults.size()];
		for (int c=0; c<energies.length ; c++)
			if (allResults.get(c).evEnergy>medianEV)
				energies[c] = allResults.get(c).score + largeEnergyAddition; 
			else
				energies[c] = allResults.get(c).score;
		int[] sortedScores = findTopMinArray(energies, takeForSCMOD, infiniteEnergy);


		// Third: Running SCMOD and saving
		for (int modelNum=0 ; modelNum<sortedScores.length ; modelNum++) {
			// Restoring the native environment
			for (int atomCounter=0 ; atomCounter<prot.atoms().size() ; atomCounter++) 
				if ((prot.atoms().atomAt(atomCounter).residueNumber()<resStart) || (prot.atoms().atomAt(atomCounter).residueNumber()>resEnd)) {
					Atom atom = protCopy.findAtomInList(prot.atoms().atomAt(atomCounter).name(), prot.atoms().atomAt(atomCounter).residueNumber());
					prot.atoms().atomAt(atomCounter).setXYZ(atom.x(), atom.y(), atom.z());
				}

			// Putting the loop in
			restoreLoopCoordinates( allResults.get(sortedScores[modelNum]).coors);
			restorePhiPsiOfLoop(allResults.get(sortedScores[modelNum]).myPhiPsi);

			// Defrosting only close residues
			prot.freeze();
			defrostCloseResidues(DEFROST_DIST);

			//"SCMODing"
			SCMOD.scmod(commands , rotLib, prot , 2, 0.0, pp);    	
			try {
				prot.atoms().print(new MeshiWriter(writePath+"/"+modelNum+".pdb"));
				if (ref!=null)
					System.out.println("888888 " + modelNum + " " + energies[sortedScores[modelNum]] + " " + calcRMS(resStart, resEnd));
				else 
					System.out.println("888888 " + modelNum + " " + energies[sortedScores[modelNum]] + " -1");
			} catch (IOException e) {
				System.out.println("No luck in writing:" + writePath+"/"+modelNum+".pdb");
			}	
		}

		System.out.print("\n\nENDENDEND");
	}
	
	

//All residues whose Cb is closer than this to the loop, 
//will be defrosted.
private void defrostCloseResidues(double defrost_dist) {
	Atom atom, loopAtom;
	Residue res;
	for (int atomCounter=0 ; atomCounter<prot.atoms().size(); atomCounter++) {
		atom = prot.atoms().atomAt(atomCounter);
		if (atom.name().equals("CB") && atom.frozen()) {
			for (int loopResidueCounter=resStart ; (loopResidueCounter<=resEnd) && (atom.frozen()); loopResidueCounter++) {
				res = prot.residue(loopResidueCounter);
				for (int loopAtomCounter=0 ; loopAtomCounter<res.atoms().size(); loopAtomCounter++) {
					loopAtom = res.atoms().atomAt(loopAtomCounter);
					if (((loopAtom.x()-atom.x())*(loopAtom.x()-atom.x()) +
							(loopAtom.y()-atom.y())*(loopAtom.y()-atom.y()) +
							(loopAtom.z()-atom.z())*(loopAtom.z()-atom.z())) < (defrost_dist*defrost_dist)) {
						atom.residue().defrost();
					}
				}
			}
		}
	}
}


//Returns the global RMS between the current loop conformation and native 
protected double calcRMS(int start, int end) {
	double totRms = 0.0;
	int counter = 0;
	for (int c=start; c<=end ; c++) {
		Atom atom = prot.residue(c).atoms().findAtomInList("N",c);
		totRms += (atom.x() - refCoors[0][3*(c-resStart)+0])*
		(atom.x() - refCoors[0][counter]) +
		(atom.y() - refCoors[1][counter])*
		(atom.y() - refCoors[1][counter]) + 
		(atom.z() - refCoors[2][counter])*
		(atom.z() - refCoors[2][counter]);
		counter++;
		atom = prot.residue(c).atoms().findAtomInList("CA",c);
		totRms += (atom.x() - refCoors[0][counter])*
		(atom.x() - refCoors[0][counter]) +
		(atom.y() - refCoors[1][counter])*
		(atom.y() - refCoors[1][counter]) + 
		(atom.z() - refCoors[2][counter])*
		(atom.z() - refCoors[2][counter]);
		counter++;
		atom = prot.residue(c).atoms().findAtomInList("C",c);
		totRms += (atom.x() - refCoors[0][counter])*
		(atom.x() - refCoors[0][counter]) +
		(atom.y() - refCoors[1][counter])*
		(atom.y() - refCoors[1][counter]) + 
		(atom.z() - refCoors[2][counter])*
		(atom.z() - refCoors[2][counter]);
		counter++;
	}
	return Math.sqrt(totRms/counter);
}

//Storing the reference loop atom coordinates in a local array 
protected void updateRefCoors() {
	refCoors = new double[3][gapLength*3];
	Atom atom;
	int counter = 0;
	for (int c=resStart ; c<=resEnd ; c++) {
		atom = ref.atoms().findAtomInList("N" , c);
		refCoors[0][counter] = atom.x();			
		refCoors[1][counter] = atom.y();			
		refCoors[2][counter] = atom.z();			
		counter++;
		atom = ref.atoms().findAtomInList("CA" , c);
		refCoors[0][counter] = atom.x();			
		refCoors[1][counter] = atom.y();			
		refCoors[2][counter] = atom.z();			
		counter++;
		atom = ref.atoms().findAtomInList("C" , c);
		refCoors[0][counter] = atom.x();			
		refCoors[1][counter] = atom.y();			
		refCoors[2][counter] = atom.z();
		counter++;
	}	
}

protected double[][] saveLoopCoordinatesCaCb() {
	double[][] savedLoopCoordinates = new double[3][gapLength*2];
	Atom atom;
	int counter=0;
	for (int c=resStart ; c<=resEnd ; c++) 
		for (int cc=0 ; cc<prot.residue(c).atoms().size() ; cc++){
			atom = prot.residue(c).atoms().atomAt(cc);
			if (atom.name().equals("CA") || atom.name().equals("CB")) {
				savedLoopCoordinates[0][counter] = atom.x();			
				savedLoopCoordinates[1][counter] = atom.y();			
				savedLoopCoordinates[2][counter] = atom.z();
				counter++;
			}
		}	
	return savedLoopCoordinates;
}

protected double[][] saveLoopCoordinatesHeavyBackbone() {
	double[][] savedLoopCoordinates = new double[3][gapLength*5];
	Atom atom;
	int counter=0;
	for (int c=resStart ; c<=resEnd ; c++) 
		for (int cc=0 ; cc<prot.residue(c).atoms().size() ; cc++){
			atom = prot.residue(c).atoms().atomAt(cc);
			if (atom.isBackbone && !atom.isHydrogen) {
				savedLoopCoordinates[0][counter] = atom.x();			
				savedLoopCoordinates[1][counter] = atom.y();			
				savedLoopCoordinates[2][counter] = atom.z();
				counter++;
			}
		}	
	return savedLoopCoordinates;
}


protected double[][] saveLoopCoordinatesNCaC() {
	double[][] savedLoopCoordinates = new double[3][gapLength*3];
	Atom atom;
	int counter = 0;
	for (int c=resStart ; c<=resEnd ; c++) {
		atom = prot.residue(c).atoms().findAtomInList("N",c);
		savedLoopCoordinates[0][counter] = atom.x();
		savedLoopCoordinates[1][counter] = atom.y();
		savedLoopCoordinates[2][counter] = atom.z();
		counter++;
		atom = prot.residue(c).atoms().findAtomInList("CA",c);
		savedLoopCoordinates[0][counter] = atom.x();
		savedLoopCoordinates[1][counter] = atom.y();
		savedLoopCoordinates[2][counter] = atom.z();
		counter++;
		atom = prot.residue(c).atoms().findAtomInList("C",c);
		savedLoopCoordinates[0][counter] = atom.x();
		savedLoopCoordinates[1][counter] = atom.y();
		savedLoopCoordinates[2][counter] = atom.z();
		counter++;
	}
	return savedLoopCoordinates;
}

protected double[][] saveLoopCoordinates() {
	double[][] savedLoopCoordinates = new double[3][gapLength*15];
	Atom atom;
	int counter=0;
	for (int c=resStart ; c<=resEnd ; c++) 
		for (int cc=0 ; cc<prot.residue(c).atoms().size() ; cc++){
			atom = prot.residue(c).atoms().atomAt(cc);
			savedLoopCoordinates[0][counter] = atom.x();			
			savedLoopCoordinates[1][counter] = atom.y();			
			savedLoopCoordinates[2][counter] = atom.z();
			counter++;
		}	
	return savedLoopCoordinates;
}

protected void restoreLoopCoordinates(double[][] savedLoopCoordinates) {
	Atom atom;
	int counter=0;
	for (int c=resStart ; c<=resEnd ; c++) 
		for (int cc=0 ; cc<prot.residue(c).atoms().size() ; cc++){
			atom = prot.residue(c).atoms().atomAt(cc);
			atom.setXYZ(savedLoopCoordinates[0][counter],savedLoopCoordinates[1][counter],
					savedLoopCoordinates[2][counter]);
			counter++;
		}
}

//This method translates the initial positions of the loop atoms away so they don't interact with the protein at all.
//Non-interaction is achieved by having no EV energy on the loop. The 12.0 Angs is just a thumb figure.
protected void putInitialLoopFar() {
	if (ref != null)
		System.out.println("The initial RMS: " + calcRMS(resStart, resEnd) );
	do {
		for (int c=resStart ; c<=resEnd ; c++) {
			for (int cc=0 ; cc<prot.residue(c).atoms().size() ; cc++) {
				prot.residue(c).atoms().atomAt(cc).addToX(12.0);
			}
		}
		try {
			energy.update();
		} 
		catch (UpdateableException e) {
			throw new RuntimeException(e);
		}	
	} while (energyTermEVloopAssembly.evaluate()>0.0001);
}

//Update the sequence array, and find which residues are pre-proline
protected void updateSeqAndPrePRO() {
	seq = Protein.getSeqOfProt(prot,resStart,resEnd);

	// Finding pre-Prolines in the sequence
	prepro = new boolean[gapLength];
	int[] tmpSEQ = Protein.getSeqOfProt(prot,resStart,resEnd+1);
	for (int c=0 ; c<gapLength ; c++)
		if ((tmpSEQ[c+1]==PRO) && (tmpSEQ[c]!=GLY) && (tmpSEQ[c]!=PRO))
			prepro[c] = true;
		else 
			prepro[c] = false;
}

protected void updatePhiPsi(int indInCorpus, int Rstart, int Rend) {
	for (int res=Rstart ; res<=Rend ; res++) {
		pp[res][0] = corpus.torsions[indInCorpus+res-Rstart][1]; // Transfering Phi
		pp[res][1] = corpus.torsions[indInCorpus+res-Rstart][2]; // Transfering Psi
/* I changed this part on 12.5.2009 after the table-propensity proved to be so good in the fragment picking.
		if (prepro[res-resStart])
			props[res-resStart] = corpus.prePro[indInCorpus+res-Rstart];
		else
			props[res-resStart] = corpus.energies[indInCorpus+res-Rstart][0][seq[res-resStart]]; 		
*/
//		if (prepro[res-resStart])
//			props[res-resStart] = corpus.propensityAll.preProVal(pp[res][0], pp[res][1]);
//		else
//			props[res-resStart] = corpus.propensityAll.propVal(seq[res-resStart], pp[res][0], pp[res][1]); 		
	}
}


protected double[][] savePhiPsiOfLoop() {
	double[][] results = new double[gapLength][2];
	for (int c=0 ; c<gapLength ; c++) {
		results[c][0] = pp[resStart+c][0];
		results[c][1] = pp[resStart+c][1];
	}
	return results;
}

private void restorePhiPsiOfLoop(double[][] ppArray) {
	for (int c=0 ; c<gapLength ; c++) {
		pp[resStart+c][0] = ppArray[c][0];
		pp[resStart+c][1] = ppArray[c][1];
	}
}

/*
 * [   ][   ][breakStart]--------[breakEnd][   ][   ]
 * 'breakStart' is the residue number of the last residue before the break.
 * 'breakEnd' is the residue number of the first residue after the break.
 */
protected double temporaryDistEnergy(int breakStart, int breakEnd) {
	double[][] data = 
	{{    0.0000,   5.9000,   6.0000,   6.1000,   6.1500},
			{    1.0000,   5.9000,   6.0000,   6.1000,   6.1500},
			{   2.0000,   5.9000,   6.0000,   6.1000,   6.1500},
			{    3.0000,   9.0000,   9.4000,   9.6000,   9.7000},
			{    4.0000,  12.0000,  12.7000,  13.0000,  13.1000},
			{    5.0000,  14.9000,  15.9000,  16.3000,  16.5000},
			{    6.0000,  17.5000,  19.0000,  19.6000,  19.9000},
			{    7.0000,  19.7000,  22.1000,  22.8000,  23.1000},
			{    8.0000,  21.7000,  24.9000,  25.8500,  26.4000},
			{    9.0000,  23.4000,  27.7000,  28.9000,  29.4000},
			{   10.0000,  24.8000,  30.1000,  31.7000,  32.4000},
			{   11.0000,  25.9000,  32.2000,  34.5000,  35.4500},
			{   12.0000,  26.8000,  34.1000,  37.1500,  38.2000},
			{   13.0000,  27.4000,  35.7500,  39.5000,  40.9000},
			{   14.0000,  28.1000,  37.1000,  41.5000,  43.3000},
			{   15.0000,  28.6000,  38.0000,  43.0500,  45.6000},
			{   16.0000,  29.0000,  39.0000,  44.5000,  48.4000},
			{   17.0000,  29.5000,  39.7500,  45.9000,  50.4500},
			{   18.0000,  29.9000,  40.4000,  47.0000,  52.8000},
			{   19.0000,  30.4000,  41.1000,  47.9000,  54.2000}};
	int len =  breakEnd-breakStart-1;
	if (len>19)
		len = 19;
	double dist = prot.residue(breakEnd).atoms().findAtomInList("N" , breakEnd).distanceFrom(
			prot.residue(breakStart).atoms().findAtomInList("C" , breakStart)) - 2.4;  // The 2.4 is the length of two peptide bonds
	if (dist>data[len][4])
		return 10.0;
	if (dist>data[len][3])
		return 4.3 + (dist-data[len][3])/(data[len][4]-data[len][3])*2.3;
	if (dist>data[len][2])
		return 2.3 + (dist-data[len][2])/(data[len][3]-data[len][2])*2.3;
	if (dist>data[len][1])
		return (dist-data[len][1])/(data[len][2]-data[len][1])*2.3;
	return 0.0;
}

/**
 * This is the recursive method, that actually build the loop. The method assume that fragments from the
 * local first (libCounter) libs has been inserted. It then continue to add fragments from the next lib, 
 * or if it is the last lib then it will analyze for closure and loop saving. When the method exit it will
 * return the coordinates of the loop to what they were before.
 * 
 * @param libCounter - see above.
 */
protected void addFromLib(int libCounter) {
	if (libCounter==0) {  // On the begining of the recursion
		loopCompleted = false;
	}
	double[][] tmpCoor = saveLoopCoordinates();
	int taken = 0;
	int libFragToChoose = 0;
	for (int c=0 ; !loopCompleted && continueSearch(libCounter,taken,c) ; c++ ) {
		libFragToChoose = chooseFrag(libCounter,c);
		if (libManners[libCounter]==-1) { // From N
			if (libs[libCounter].isFragCompatible(libFragToChoose,resStart+libStarts[libCounter]-libs[libCounter].overlap(),libs[libCounter].overlap(),-1,similarityMatchCO)) {
				libs[libCounter].insertFragToProt(rotLib, libFragToChoose, prot, resStart+libStarts[libCounter]-libs[libCounter].overlap(), libs[libCounter].overlap(), -1);
				updatePhiPsi(libs[libCounter].libOrig(libFragToChoose)+libs[libCounter].overlap(), resStart+libStarts[libCounter], resStart+libEnds[libCounter]);
				fragRank[libCounter] = libFragToChoose;
				if (libCounter==(libs.length-1)) { // closing
					if (loopClosingCondition()) { 						
						taken++;
						actuallyTaken[libCounter] = taken;
						analyzedCloseLoop();
						loopCompleted = true;
					}
				}
				else {
					if (continueForNextCall(libCounter)) {
						taken++;
						actuallyTaken[libCounter] = taken;
						addFromLib(libCounter+1);
					}
				}
			}
		}
		else { // From C
			if (libs[libCounter].isFragCompatible(libFragToChoose,resStart+libStarts[libCounter],libs[libCounter].overlap(),1,similarityMatchCO)) {
				libs[libCounter].insertFragToProt(rotLib, libFragToChoose, prot, resStart+libStarts[libCounter], libs[libCounter].overlap(), 1);
				updatePhiPsi(libs[libCounter].libOrig(libFragToChoose), resStart+libStarts[libCounter], resStart+libEnds[libCounter]);
				fragRank[libCounter] = libFragToChoose;
				if (libCounter==(libs.length-1)) { // closing
					if (loopClosingCondition()) { 						
						taken++;
						actuallyTaken[libCounter] = taken;
						analyzedCloseLoop();
						loopCompleted = true;
					}
				}
				else {
					if (continueForNextCall(libCounter)) {
						taken++;
						actuallyTaken[libCounter] = taken;
						addFromLib(libCounter+1);
					}
				}
			}					
		}
	}
	restoreLoopCoordinates(tmpCoor);
	return;
}

protected void resetTotalEnergy() {
	// For the sake of calculating HB fast, only the loop is activated.
	for (int c=0; c<prot.atoms().size() ; c++)
		if ((prot.atoms().atomAt(c).residueNumber()>=resStart) && 
				(prot.atoms().atomAt(c).residueNumber()<=resEnd)) 
			prot.atoms().atomAt(c).activate();
		else 
			prot.atoms().atomAt(c).inActivate();
	// Constructing the total energy 
	EnergyCreator[] energyCreators = getEnergyCreatorArray();
	energy = new TotalEnergy(prot, new DistanceMatrix(prot.atoms(),  5.5, 2.0, 4), energyCreators, commands);
	energyTermHB = (SimpleHydrogenBondEnergy) (energy.getEnergyTerm(new SimpleHydrogenBondEnergy()));
	energyTermEVloopFinal = (SoftExcludedVol) (energy.getEnergyTerms(new SoftExcludedVol())[0]);
	energyTermEVloopAssembly = (SoftExcludedVol) (energy.getEnergyTerms(new SoftExcludedVol())[1]);
	((SoftExcludedVolEnergyElement) (energyTermEVloopFinal.getElement())).setFragLimits(resStart , resEnd);	
	((SoftExcludedVolEnergyElement) (energyTermEVloopAssembly.getElement())).setFragLimits(resStart , resEnd);
	try {
		energy.update();
	} 
	catch (UpdateableException e) {
		throw new RuntimeException(e);
	}	
	System.out.println("Initial EV energy of the loop is: " + energyTermEVloopFinal.evaluate() + "\n" +
			"Initial HB energy of the loop is: " + energyTermHB.evaluate());
	setUpLocalEV();
}

protected EnergyCreator[] getEnergyCreatorArray() {
	EnergyCreator[] energyCreators = {new  SoftExcludedVolCreator(1.0 , 9 , 1.0),
			new  SoftExcludedVolCreator(1.0 , 10 , 1.0),
			new SimpleHydrogenBond_Dahiyat_LowAccuracy_BBonly_Creator(1.0)};
	return energyCreators;
}


protected double evaluateEVextLoop() {
	return evaluateEV(evData_extLoop);
}

protected double evaluateEVinterLoop() {
	return evaluateEV(evData_interLoop);
}

protected double evaluateTouchLoop() {
	return evaluateTouch(touchData_extLoop);
}

private double evaluateEV(Vector<LocalEVdata> vec) {
	double ev = 0;
	double tmpRsquared;
	double tmpR;
	for (LocalEVdata evData : vec) {
		tmpRsquared = (evData.atom1.x()-evData.atom2.x())*
		              (evData.atom1.x()-evData.atom2.x()) + 
		              (evData.atom1.y()-evData.atom2.y())*
		              (evData.atom1.y()-evData.atom2.y()) +
		              (evData.atom1.z()-evData.atom2.z())*
		              (evData.atom1.z()-evData.atom2.z());
		if (tmpRsquared<evData.rLJsquared) {
			tmpR = Math.sqrt(tmpRsquared);
			ev += (tmpR-evData.rLJ)*(tmpR-evData.rLJ);
		}
	}
	return ev;
}

private double evaluateTouch(Vector<LocalEVdata> vec) {
	double touch = 0;
	double tmpRsquared;
	for (LocalEVdata evData : vec) {
		tmpRsquared = (evData.atom1.x()-evData.atom2.x())*
		              (evData.atom1.x()-evData.atom2.x()) + 
		              (evData.atom1.y()-evData.atom2.y())*
		              (evData.atom1.y()-evData.atom2.y()) +
		              (evData.atom1.z()-evData.atom2.z())*
		              (evData.atom1.z()-evData.atom2.z());
		if (tmpRsquared<evData.rLJsquared) {
			touch++;
		}
	}
	return touch;
}

private void setUpLocalEV() {
	SoftExcludedVolParametersList parametersList = 
		new SoftExcludedVolParametersList(commands.firstWord(PARAMETERS_DIRECTORY).secondWord() + "/" + 
				EXCLUDED_VOL_PARAMETERS,1.0);	
	// Calculating the Inter-loop data
	Atom atom1,atom2;
	double rLJ;
	evData_interLoop = new Vector<LocalEVdata>();
	for (int c1=0 ; c1<prot.atoms().size() ; c1++)
		for (int c2=c1+1 ; c2<prot.atoms().size() ; c2++) {
			atom1 = prot.atoms().atomAt(c1);
			atom2 = prot.atoms().atomAt(c2);
			if (!atom1.isHydrogen && !atom2.isHydrogen && 
					atom1.isBackbone && atom2.isBackbone &&
					atomInLoop(atom1) && atomInLoop(atom2) &&
					!inBondedList(atom1, atom2)) {
				rLJ = ((SoftExcludedVolParameters) parametersList.parameters(new Distance(atom1,atom2))).sigma;
				evData_interLoop.add(new LocalEVdata(atom1,atom2,rLJ));
			}
		}
	// Calculating the ext-loop data
	evData_extLoop = new Vector<LocalEVdata>();
	for (int c1=0 ; c1<prot.atoms().size() ; c1++)
		for (int c2=c1+1 ; c2<prot.atoms().size() ; c2++) {
			atom1 = prot.atoms().atomAt(c1);
			atom2 = prot.atoms().atomAt(c2);
			if (!atom1.isHydrogen && !atom2.isHydrogen && 
					atom1.isBackbone && atom2.isBackbone &&
					((!atomInLoop(atom1) && atomInLoop(atom2)) ||
					(atomInLoop(atom1) && !atomInLoop(atom2))) &&
					!inBondedList(atom1, atom2)){
				rLJ = ((SoftExcludedVolParameters) parametersList.parameters(new Distance(atom1,atom2))).sigma;
				evData_extLoop.add(new LocalEVdata(atom1,atom2,rLJ));
			}
		}	
	// Calculating the data for touches of the loop with the outer loop
	touchData_extLoop = new Vector<LocalEVdata>();
	for (int c1=0 ; c1<prot.atoms().size() ; c1++)
		for (int c2=c1+1 ; c2<prot.atoms().size() ; c2++) {
			atom1 = prot.atoms().atomAt(c1);
			atom2 = prot.atoms().atomAt(c2);
			if (atom1.name().equals("CB") && atom2.name().equals("CB") &&
					((!atomInLoop(atom1) && atomInLoop(atom2)) ||
					(atomInLoop(atom1) && !atomInLoop(atom2))) &&
					!inBondedList(atom1, atom2)){
				touchData_extLoop.add(new LocalEVdata(atom1,atom2,  6.0  ));
			}
		}	

}

protected boolean atomInLoop(Atom atom) {
	return ((atom.residueNumber()>=resStart) && (atom.residueNumber()<=resEnd));
}

protected boolean inBondedList(Atom atom1 , Atom atom2) {
	return energy.distanceMatrix().bondedList().contains(new Distance(atom1,atom2));
}


//Here we take at most N (could be less) indices whose values are the lowest AND less than the maximalVal 
public static int[] findTopMinArray(double[] ar, int N, double maximalVal) {
	int[] result = new int[N]; 
	boolean[] notVisit = new boolean[ar.length];
	for (int c=0; c<ar.length ; c++) 
		notVisit[c] = true;
	int n;
	for (n=0 ; n<N ; n++) {
		double minVal = 1e308;
		int minValInd = -1;
		for (int m=0 ; m<ar.length ; m++)
			if (notVisit[m] && (ar[m]<minVal) && (ar[m]<maximalVal)) {
				minVal = ar[m];
				minValInd = m;
			}
		if (minValInd>-1) {
			result[n] = minValInd;
			notVisit[minValInd] = false;
		}
		else
			break;
	}
	if (n==N)
		return result;
	else {
		int[] shorterResults = new int[n];
		for (int nn=0; nn<n; nn++)
			shorterResults[nn] = result[nn];
		return shorterResults;
	}
}


protected abstract void setClosureAtoms();

protected abstract void buildCoarsFragments();

//Build the local frag libraries. 
protected abstract void buildLibraries(); 

/**
 * This method will give the suitable fragment length to work with.
 */
protected abstract int giveFragLength();

/**
 * This method is called from addToLib. It should decide whether to continue the search from this lib.
 * @param libCounter - The current local lib being searched
 * @param taken - how many were already analyzed from this lib.
 * @param c - the local counter in the call of addToLib
 */
protected abstract boolean continueSearch(int libCounter, int taken, int c);

/**
 * This method chose a fragment from the local frag library to be tested in the loop.
 * @param libCounter - The current local lib being searched
 * @param c - the local counter in the call of addToLib
 */
protected abstract int chooseFrag(int libCounter, int c);

/**
 * This method decides if we can recursively continue, after the last fragment was inserted.
 * @param libCounter - The current local lib being searched
 */
protected abstract boolean continueForNextCall(int libCounter);

/**
 * This method checks if some condition upon loop closing is fullfilled. If so, it will continue to 
 * analyze the loop.
 */
protected abstract boolean loopClosingCondition();

/**
 * What to do when the loop is formed.
 */
protected abstract void analyzedCloseLoop();


// This class is for the local EV evaluation
private class LocalEVdata {
	protected Atom atom1;
	protected Atom atom2;
	protected double rLJ;
	protected double rLJsquared;
	
	public LocalEVdata(Atom atom1, Atom atom2, double rLJ) {
		this.atom1 = atom1;
		this.atom2 = atom2;
		this.rLJ = rLJ;
		rLJsquared = rLJ*rLJ;
	}
}


} // Of Loop



