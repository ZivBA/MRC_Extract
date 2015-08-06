package meshi.applications.loopBuilding.applications;

import java.util.Calendar;
import java.util.StringTokenizer;
import java.util.Vector;

import meshi.applications.corpus.Corpus;
import meshi.applications.corpus.PropensityMatrix;
import meshi.applications.loopBuilding.LoopBuilderDevelopePhase2SecondAttempt;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.softExcludedVol.SoftExcludedVolCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LBFGS;
import meshi.optimizers.Minimizer;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.MeshiProgram;
import meshi.util.Isite.IsiteLib;
import meshi.util.file.File2StringArray;
import programs.PutHydrogens;

/**
 * 
 * @author Nir
 * 
 * This class will let you chose the frags to model.
 * The numbering of the type are:
 * -1,1 - regular N- C- term frags, connecting to protein.
 * -2,2 - like above but NOT connecting to a known structure.
 * -3,3 - like above but fixed.
 * <-999, >999 , like above but interfacing a fixed frags on the other side. The frag index is found by (X-1000) or (-1000-X) depending on the sign. 
 * 
 */
public class DoBuildCoarseLoopWithClosure extends MeshiProgram implements Residues, 
AtomTypes{ /**
 * The implemented
 * interfaces defines the 
 * names of atom and residue 
 * types. 
 **/

	/**
	 * A string with the name of the pdb file to minimize.
	 **/

	private static CommandList commands; 

	private static String fileNameCorpus = null;

	private static String fileNameProt = null;

	private static String fileNameRef = null;

	private static String fileNameHomologues = null;

	private static String pathString = null;

	private static int resStart = -1;

	private static int resEnd = -1;

	private static int libSize = -1;

	private static double rmsCutOff = -1.0;

	private static double rmsMatchCO = -1.0;

	private static double closureTolerance = -1.0;

	private static double TAKE_AROUND_LOOP = 10.0;    

	private static Vector<int[]> frags; 

	private static int maxNumOfCoarse = Integer.MAX_VALUE; 

	private static int numToPrint = 2500; 

	public static void main(String[] args) {
		init(args);
		// Building the corpus
		Corpus corpus = new Corpus(fileNameCorpus);
//		Use this when you want known fragments in your fragment library
//		---------------------------------------------------------------
//		EnergyCreator[] tmpEnergyCreators = { };
//		CorpusQuickNoSolvRot1 tmpCorpus = new CorpusQuickNoSolvRot1(fileNameProt, commands,
//				tmpEnergyCreators);
//		corpus.merge(tmpCorpus);
//		---------------------------------------------------------------
		String[] excludedProts = getExcludedProt(fileNameHomologues);
		corpus.setExcludedProteinsFromUngapped(excludedProts);
		corpus.propensityAll = new PropensityMatrix(commands,"allResPropensity.dat");
		corpus.propensityLoop = new PropensityMatrix(commands,"loopResPropensity.dat");

		
		// Highly truncated proteins for loop building:
		// --------------------------------------------	
		// Building the protein
		Protein tmpprot = new Protein(new AtomList(fileNameProt) , new ResidueExtendedAtoms(ADD_ATOMS));
		PutHydrogens.adjustBackboneStereochemistry(commands, tmpprot, resStart, resEnd);
		PutHydrogens.adjustHydrogens(commands, tmpprot);
		Protein prot = takeNearResidues(tmpprot,resStart, resEnd, 12.0);
		for (int cc=0 ; cc<prot.atoms().size() ; cc++)
			prot.atoms().atomAt(cc).setChain("A");
		// Building the Reference
		tmpprot = new Protein(new AtomList(fileNameRef) , new ResidueExtendedAtoms(ADD_ATOMS));
		PutHydrogens.adjustHydrogens(commands, tmpprot);
		Protein ref = takeNearResidues(tmpprot,resStart, resEnd, 13.0);
		for (int cc=0 ; cc<ref.atoms().size() ; cc++)
			ref.atoms().atomAt(cc).setChain("A");
		
		// Building the Isite Library
		String parametersDirectory = commands.firstWord(KeyWords.PARAMETERS_DIRECTORY).secondWord();
		IsiteLib iLib = new IsiteLib(parametersDirectory +"/Isite51");
		iLib.setCorpus(corpus);
		iLib.loadMotifInCorpus(parametersDirectory +"/motifAssignment.dat");
		iLib.loadAndAnalyzeSequenceIntoLib(Protein.getSeqOfProt(tmpprot, 0, tmpprot.lastResidue()));
		iLib.loadAndAnalyzeBackboneIntoLib(tmpprot);
		for (int pos=0 ; pos<tmpprot.lastResidue() ; pos++) {
			iLib.bestMotifInSeqPos(pos);
		}
		corpus.iSiteLib = iLib;

		
		// weakly truncated proteins for loop building:
		// --------------------------------------------	
		// Building the protein
		tmpprot = new Protein(new AtomList(fileNameProt) , new ResidueExtendedAtoms(ADD_ATOMS));
		PutHydrogens.adjustHydrogens(commands, tmpprot);
		Protein evalProt = takeNearResidues(tmpprot,resStart, resEnd, TAKE_AROUND_LOOP);
		for (int cc=0 ; cc<evalProt.atoms().size() ; cc++)
			evalProt.atoms().atomAt(cc).setChain("A");		

		// "Looping"
		makeFragsDescription();
		long startingTime = Calendar.getInstance().getTimeInMillis();
		
		//	For Phase 1 in the second round of tries:
		//  -------------
//		LoopBuilderDevelopePhase1SecondAttempt loop = new LoopBuilderDevelopePhase1SecondAttempt(
//				commands, pathString, corpus, prot, ref,
//				resStart, resEnd, rmsMatchCO, rmsCutOff,
//				frags, closureTolerance); 
		//	For Phase 2 in the second round of tries:
		//  -------------
		LoopBuilderDevelopePhase2SecondAttempt loop = new LoopBuilderDevelopePhase2SecondAttempt(
				commands, pathString, corpus, prot, ref,
				resStart, resEnd, rmsMatchCO, rmsCutOff,
				frags, closureTolerance); 
		loop.runCoarse(maxNumOfCoarse,numToPrint);
		long endTime = Calendar.getInstance().getTimeInMillis();
		System.out.println("Time for creating loopsoop: " + ((endTime-startingTime)/1000) + " seconds.");
		System.out.println("12121212 " + ((endTime-startingTime)/1000));
		loop.setProt(evalProt);
//		loop.fragmentPickingAnalysis(0);
		loop.developementAnalysis(0);
		loop.setProt(prot);
		endTime = Calendar.getInstance().getTimeInMillis();
		System.out.println("Net runtime of loop: " + ((endTime-startingTime)/1000) + " seconds.");
		System.out.println("13131313 " + ((endTime-startingTime)/1000));

		//	For Phase 1 in the first round of tries (when I was trying to do the ROT1 potential):
		//  -------------
		//		LoopBuilderDevelopePhase1OpenCoarseLoops loop = new LoopBuilderDevelopePhase1OpenCoarseLoops(
		//				commands, pathString, corpus, prot, ref,
		//				resStart, resEnd, rmsMatchCO, rmsCutOff,
		//				frags, closureTolerance); 
		//		for (int module = 0 ; module<maxNumOfCoarse ; module+=1000) {
		//			loop.runCoarse(999,numToPrint);
		//			loop.setProt(evalProt);
		//			loop.developementAnalysis(module);
		//			loop.setProt(prot);
		//		}

		//	For Phase 2:
		//  -------------
		//		LoopBuilderDevelopePhase2OpenCoarseLoops loop = new LoopBuilderDevelopePhase2OpenCoarseLoops(
		//				commands, pathString, corpus, prot, ref,
		//				resStart, resEnd, rmsMatchCO, rmsCutOff,
		//				frags, closureTolerance); 
		//		int runningNumber = 0;
		//		for (int module = 0 ; module<maxNumOfCoarse ; module+=1000) {
		//			loop.runCoarse(999,numToPrint);
		//			loop.developementAnalysis(runningNumber);
		//			runningNumber += loop.getAllResults().size();
		//		}
	}

	private static String[] getExcludedProt(String fileName) {
		String[] tmp = File2StringArray.f2a(fileName);
		String[] results = new String[tmp.length];
		for (int c=0 ; c<tmp.length ; c++) {
			StringTokenizer st = new StringTokenizer(tmp[c]);
			results[c] = st.nextToken().trim();
			System.out.println("Excluding from corpus protein: " + results[c]);
		}
		return results;
	}

	private static void makeFragsDescription() {
		// ****************************************
		// In the first try the DEFAULT_SIZE was 3.
		// ****************************************
		int DEFAULT_SIZE = 2;
		frags = new Vector<int[]>();
		// Two anchors
		int[] anchorN = {resStart, resStart+DEFAULT_SIZE-1, -1, libSize/2};
		int[] anchorC = {resEnd-DEFAULT_SIZE+1, resEnd, 1, libSize/2};
		System.out.println("Fragment 1: " + anchorN[0] + "-" + anchorN[1] + " lib size:" + anchorN[3] + "    type:" + anchorN[2]);
		System.out.println("Fragment 2: " + anchorC[0] + "-" + anchorC[1] + " lib size:" + anchorC[3] + "    type:" + anchorC[2]);
		frags.add(anchorN);
		frags.add(anchorC);

		int fragCounter = 3;
		int startRes = resStart+DEFAULT_SIZE;
		int endRes = resEnd-DEFAULT_SIZE;
		int manner = -2;
		while (startRes<=endRes) {
			int fragSize = Math.min(DEFAULT_SIZE, endRes-startRes+1);
			if (manner==-2) {
				int[] tmp = {startRes, startRes+fragSize-1, -2, libSize};
				frags.add(tmp);
				startRes += fragSize;
			}
			else {
				int[] tmp = {endRes-fragSize+1, endRes, 2, libSize};
				frags.add(tmp);
				endRes -= fragSize;				
			}
			System.out.println("Fragment " + fragCounter + ": " + frags.lastElement()[0] + "-" + frags.lastElement()[1] + " lib size:" + frags.lastElement()[3] + "    type:" + frags.lastElement()[2]);
			fragCounter++;
			manner = -manner;
		}
	}

	protected static Protein takeNearResidues(Protein tmpprot, int start, int end, double disCO) {
		tmpprot.freeze();
		for (int c=start; c<=end ; c++)
			tmpprot.residue(c).atoms().defrost();
		for (int c=0; c<tmpprot.atoms().size() ; c++) 
			for (int d=0; d<tmpprot.atoms().size() ; d++) 
				if (((tmpprot.atoms().atomAt(c).x()-tmpprot.atoms().atomAt(d).x())*
						(tmpprot.atoms().atomAt(c).x()-tmpprot.atoms().atomAt(d).x()) +
						(tmpprot.atoms().atomAt(c).y()-tmpprot.atoms().atomAt(d).y())*
						(tmpprot.atoms().atomAt(c).y()-tmpprot.atoms().atomAt(d).y()) +
						(tmpprot.atoms().atomAt(c).z()-tmpprot.atoms().atomAt(d).z())*
						(tmpprot.atoms().atomAt(c).z()-tmpprot.atoms().atomAt(d).z())) < (disCO*disCO))
					if ((tmpprot.atoms().atomAt(d).residueNumber()>=start) && 
							(tmpprot.atoms().atomAt(d).residueNumber()<=end))
						tmpprot.residue(tmpprot.atoms().atomAt(c).residueNumber()).atoms().defrost();
		/* ***************************** 
		 * This is for the fragment pick part
		 ***************************** */
		for (int c=(start-9); c<=(end+9) ; c++)
			tmpprot.residue(c).atoms().defrost();
		return (new Protein(tmpprot.atoms().filter(new AtomList.IsDefrostedFilter()).duplicate() , new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS)));
	}
	
    public static void adjustBackboneStereochemistry(CommandList commandsParam, Protein prot, int fromRes, int toRes) {
    	prot.freeze();
    	for (int cc=fromRes ; cc<=toRes ; cc++) {
    		prot.residue(cc).defrost();
    	}
   		DistanceMatrix distanceMatrix = new DistanceMatrix(prot.atoms(),   5.5, 2.0, 4);  
   		EnergyCreator[] energyCreators = {  
    		    new BondCreator(),
    		    new AngleCreator(),
    		    new PlaneCreator(),
    		    new OutOfPlaneCreator(),
    		    new SoftExcludedVolCreator(10.0 , 12 , 1.0)
    	};
   		TotalEnergy energy = new TotalEnergy(prot, distanceMatrix, energyCreators, commandsParam);

   		Minimizer minimizer = new LBFGS(energy, 0.05 , 10000 , 100);  
    	try {
    	  System.out.println(minimizer.minimize());
    	}
    	catch (Exception e) {
    		System.out.println("\n\n\nThere was a problem in the minimization.\n");
    	}  
    }


	/** ================================= init =========================================
	 *
	 *A static function for parsing of the command line arguments and assigning the 
	 *variables commandsFileName, modelFileName and randomNumberSeed with the right inputs. Note that this
	 *static method is using parsing functions such as getOrderedArguments that are defined in MeshiProgram
	 *that MinimizeProtein inherits.
	 **/

	protected static void init(String[] args) {

		/**** NOTE *** the next two lines. Because of a BUG in the Java VM, the 
		 * interfaces "Residues" and "AtomTypes" are not loaded automatically when MinimizeProtein initialize. 
		 * For this purpose these two lines are crucial wherever these two interfaces are implemented. The user might 
		 * rightfully feel that these two lines are "black magic" programming, but happily to our knowledge this is 
		 * the only bizarre phenomenon we are aware of in meshi.
		 **/
		int zvl = ALA; // force the reading of "meshi.parameters.Residues"
		zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"


		String errorMessage = ("\n                  ******************\n"+
				"Usage java -Xmx600m DoLoop <command file> <corpora file> <prot corpus> <prot name> <resStart> <resEnd> \n"+
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));

		String commandsFileName = getOrderedArgument(args);
		if (commandsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# commandsFileName = "+commandsFileName);

		commands = new CommandList(commandsFileName);

		fileNameCorpus = getOrderedArgument(args);
		if (fileNameCorpus == null) throw new RuntimeException(errorMessage);
		System.out.println("# Corpus file name is "+fileNameCorpus);

		fileNameProt = getOrderedArgument(args);
		if (fileNameProt == null) throw new RuntimeException(errorMessage);
		System.out.println("# Prot file name is "+fileNameProt);

		fileNameRef = getOrderedArgument(args);
		if (fileNameRef == null) throw new RuntimeException(errorMessage);
		System.out.println("# Ref file name is "+fileNameRef);

		fileNameHomologues = getOrderedArgument(args);
		if (fileNameHomologues == null) throw new RuntimeException(errorMessage);
		System.out.println("# Homologue file name is "+fileNameHomologues);

		String tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		resStart = (new Integer(tmp)).intValue();
		System.out.println("# Start resisue  is "+resStart);

		tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		resEnd = (new Integer(tmp)).intValue();
		System.out.println("# end resisue  is "+resEnd);

		tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		libSize = (new Integer(tmp)).intValue();
		System.out.println("# Lib size is "+libSize);

		tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		rmsCutOff = (new Double(tmp)).doubleValue();
		System.out.println("# rmsCutOff for library building is "+rmsCutOff);

		tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		rmsMatchCO = (new Double(tmp)).doubleValue();
		System.out.println("# rmsMatchCO is "+rmsMatchCO);

		tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		closureTolerance = (new Double(tmp)).doubleValue();
		System.out.println("# Closure tolerance is "+ closureTolerance);

		tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		TAKE_AROUND_LOOP = (new Double(tmp)).doubleValue();
		System.out.println("# TAKE_AROUND_LOOP is "+TAKE_AROUND_LOOP);

		tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		maxNumOfCoarse = (new Integer(tmp)).intValue();
		System.out.println("# Maximal number of coarse loops: "+maxNumOfCoarse);

		tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		numToPrint = (new Integer(tmp)).intValue();
		System.out.println("# Number of coarse loops to write to disk: "+numToPrint);

		pathString = getOrderedArgument(args);
		if (pathString == null) throw new RuntimeException(errorMessage);
		System.out.println("# Path for output is "+pathString);

		tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		initRandom((new Integer(tmp)).intValue());
		System.out.println("# Random Seed is "+ (new Integer(tmp)).intValue());

	}
}
