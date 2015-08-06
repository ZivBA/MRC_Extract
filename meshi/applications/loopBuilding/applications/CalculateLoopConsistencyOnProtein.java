package meshi.applications.loopBuilding.applications;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;

import meshi.applications.prediction.GDTcalculator;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;
import meshi.util.file.MeshiWriter;
import meshi.util.rotamericTools.RotamericTools;

public class CalculateLoopConsistencyOnProtein extends MeshiProgram implements Residues, AtomTypes {

	private static String loopDirString = null;  
	private static String loopStartString = null;  
	private static String loopEndString = null;  

	
//	public static void main(String[] args) {
//		init(args); 
//		
//		// Paramteres:
//		// -----------
//		int loopLength = 12;
//		int scanShift = 3;
//		int takeToModel = 5;
//		int[] weights = {1,1,1,2,2,2,2,2,2,1,1,1};
//		double TH = 3.0;
//		
//		Protein template = new Protein(loopDir + "/homology.pdb", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
//		int[] resCounts = new int[template.residues().size() + 1]; 
//		int[][] gaps = MakeLimitsForScan.getScanLimits(template, loopLength, scanShift);
//		for (int gapC=0 ; gapC<gaps.length ; gapC++) {
//			String loopSubDir = loopDir + "/Loop_" + gaps[gapC][0] + "_" + gaps[gapC][1];
//			double[] bestLoop = bestRMSinLoopGroup(template, loopSubDir, takeToModel);
//			int bestModel = ((int) bestLoop[0]);
//			if (bestModel>=0) {
//				System.out.println("Using model: " + loopSubDir+"/"+bestModel+".pdb.min     RMS: " + bestLoop[1]);
//				Protein loop = new Protein(loopSubDir+"/"+bestModel+".pdb.min", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
//				updateCount(resCounts, loop, template, TH, weights);
//			}
//		}
//		printResults(resCounts,template);		
//	}

	// A version were all the template is considered!
//	public static void main(String[] args) {
//		init(args); 
//		
//		// Paramteres:
//		// -----------
//		int loopLength = 12;
//		int scanShift = 1;
//		int takeToModel = 5;
//
//		Protein template = new Protein(loopDir + "/homology.pdb", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
//		Protein conc = new Protein(loopDir + "/conc.pdb", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
//		Protein reference = new Protein(loopDir + "/native.pdb", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
//		GDTcalculator.alignBySubset(template.atoms(), conc.atoms(), 1.5);
//		GDTcalculator.alignBySubset(template.atoms(), reference.atoms(), 1.5);
//		try {
//			reference.atoms().print(new MeshiWriter(loopDir + "/nativeAligned.pdb"));
//		} catch (IOException e) {
//			throw new RuntimeException("Error writing.");
//		}
//		try {
//			conc.atoms().print(new MeshiWriter(loopDir + "/concAligned.pdb"));
//		} catch (IOException e) {
//			throw new RuntimeException("Error writing.");
//		}
//		DecimalFormat fmt2 = new DecimalFormat("0.##");
//		int[][] gaps = MakeLimitsForScan.getScanLimits(template, loopLength, scanShift);
//		for (int gapC=0 ; gapC<gaps.length ; gapC++) {
//			if ((gaps[gapC][1] - gaps[gapC][0] + 1) == loopLength) {
//				try {
//					String loopSubDir = loopDir + "/Loop_" + gaps[gapC][0] + "_" + gaps[gapC][1];
//					Protein shortTemplate = extractLoopFromTemplate(template , gaps[gapC][0] , gaps[gapC][1]);
//					Protein shortConc = extractLoopFromTemplate(conc , gaps[gapC][0] , gaps[gapC][1]);
//					double[] bestLoop = bestRMSinLoopGroup(template, loopSubDir, takeToModel);
//					int bestModel = ((int) bestLoop[0]);
//					if (bestModel>=0) {
//						Protein loop = new Protein(loopSubDir+"/0.pdb.min", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
//						double rmsNo1 = calcRMSonHonigBackbone(loop, reference);
//						double rmsTemp = calcRMSonHonigBackbone(shortTemplate, reference);
//						double rmsConc = calcRMSonHonigBackbone(shortConc, reference);
//						System.out.print("99999 " + gaps[gapC][0] + " " + gaps[gapC][1] + " " +
//								bestModel + " " + fmt2.format(bestLoop[1]) + " " + 
//								fmt2.format(rmsTemp) + " " + fmt2.format(rmsConc) + " " +
//								fmt2.format(rmsNo1) + " ");
//						double[] bestLoopToNat = bestRMSinLoopGroup(reference, loopSubDir, takeToModel);
//						int bestModelToNat = ((int) bestLoopToNat[0]);
//						System.out.println(fmt2.format(bestLoopToNat[1]) + " " + bestModelToNat);
//					}
//				}
//				catch (Exception e) {
//					// Do nothing
//				}
//			}
//		}
//	}

	// A version were only the part of the template around the loop is considered!
	public static void main(String[] args) throws Exception {
		init(args); 
		if (loopEndString!=null) {
			calcConsistencyOneLoop(loopDirString,Integer.parseInt(loopStartString),Integer.parseInt(loopEndString));
		}
		else {
			calcConsistencyAllGaps(loopDirString);
		}
		
		
	}
	
	public static void calcConsistencyOneLoop(String loopDir,int startRes, int endRes) throws Exception {
		// Paramteres:
		// -----------
		int takeToModel = 5;

		Protein origTemplate = new Protein(loopDir + "/homology.pdb", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		Protein conc = new Protein(loopDir + "/conc.pdb", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		Protein reference = new Protein(loopDir + "/native.pdb", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		DecimalFormat fmt2 = new DecimalFormat("0.##");
		String Soutput = "Error in run.";
		Protein template = takeNearResidues(origTemplate, startRes, endRes, 12.0);
		GDTcalculator.alignBySubset(template.atoms(), conc.atoms(), 1.5);
		GDTcalculator.alignBySubset(template.atoms(), reference.atoms(), 1.5);
		try {
			reference.atoms().print(new MeshiWriter(loopDir + "/nativeAligned.pdb"));
		} catch (IOException e) {
			throw new RuntimeException("Error writing.");
		}
		try {
			conc.atoms().print(new MeshiWriter(loopDir + "/concAligned.pdb"));
		} catch (IOException e) {
			throw new RuntimeException("Error writing.");
		}
		String loopSubDir = loopDir + "/Loop_" + startRes + "_" + endRes;
		Protein shortTemplate = extractLoopFromTemplate(template , startRes, endRes);
		Protein shortConc = extractLoopFromTemplate(conc , startRes, endRes);
		double[] bestLoop = bestRMSinLoopGroup(template, loopSubDir, takeToModel);
		int bestModel = ((int) bestLoop[0]);
		if (bestModel>=0) {
			Protein loop = new Protein(loopSubDir+"/0.pdb.min", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
			Protein loopProt = mergeTemplateAndLoop(template,loop);
			try {
				loopProt.atoms().print(new MeshiWriter(loopDir + "/loop0Aligned.pdb"));
			} catch (IOException e) {
				throw new RuntimeException("Error writing.");
			}
			double rmsNo1 = calcRMSonHonigBackbone(loop, reference);
			double rmsTemp = calcRMSonHonigBackbone(shortTemplate, reference);
			double rmsConc = calcRMSonHonigBackbone(shortConc, reference);
			Soutput = "99999 " + startRes + " " + endRes + " " +
			bestModel + " " + fmt2.format(bestLoop[1]) + " " + 
			fmt2.format(rmsTemp) + " " + fmt2.format(rmsConc) + " " +
			fmt2.format(rmsNo1) + " ";
			double[] bestLoopToNat = bestRMSinLoopGroup(reference, loopSubDir, takeToModel);
			int bestModelToNat = ((int) bestLoopToNat[0]);
			Soutput += fmt2.format(bestLoopToNat[1]) + " " + bestModelToNat;
		}
		System.out.println(Soutput);
	}

	
	public static void calcConsistencyAllGaps(String loopDir) throws Exception {
		// Paramteres:
		// -----------
		int loopLength = 12;
		int scanShift = 1;
		int takeToModel = 5;

		Protein origTemplate = new Protein(loopDir + "/homology.pdb", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		Protein conc = new Protein(loopDir + "/conc.pdb", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		Protein reference = new Protein(loopDir + "/native.pdb", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		DecimalFormat fmt2 = new DecimalFormat("0.##");
		BufferedWriter bw = new BufferedWriter(new FileWriter("output"));		
		int[][] gaps = MakeLimitsForScan.getScanLimits(origTemplate, loopLength, scanShift);
		for (int gapC=0 ; gapC<gaps.length ; gapC++) {
			if ((gaps[gapC][1] - gaps[gapC][0] + 1) == loopLength) {
				try {
					Protein template = takeNearResidues(origTemplate, gaps[gapC][0], gaps[gapC][1], 12.0);
					GDTcalculator.alignBySubset(template.atoms(), conc.atoms(), 1.5);
					GDTcalculator.alignBySubset(template.atoms(), reference.atoms(), 1.5);
					String loopSubDir = loopDir + "/Loop_" + gaps[gapC][0] + "_" + gaps[gapC][1];
					Protein shortTemplate = extractLoopFromTemplate(template , gaps[gapC][0] , gaps[gapC][1]);
					Protein shortConc = extractLoopFromTemplate(conc , gaps[gapC][0] , gaps[gapC][1]);
					double[] bestLoop = bestRMSinLoopGroup(template, loopSubDir, takeToModel);
					int bestModel = ((int) bestLoop[0]);
					if (bestModel>=0) {
						Protein loop = new Protein(loopSubDir+"/0.pdb.min", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
						double rmsNo1 = calcRMSonHonigBackbone(loop, reference);
						double rmsTemp = calcRMSonHonigBackbone(shortTemplate, reference);
						double rmsConc = calcRMSonHonigBackbone(shortConc, reference);
						bw.write("99999 " + gaps[gapC][0] + " " + gaps[gapC][1] + " " +
								bestModel + " " + fmt2.format(bestLoop[1]) + " " + 
								fmt2.format(rmsTemp) + " " + fmt2.format(rmsConc) + " " +
								fmt2.format(rmsNo1) + " ");
						System.out.print("99999 " + gaps[gapC][0] + " " + gaps[gapC][1] + " " +
								bestModel + " " + fmt2.format(bestLoop[1]) + " " + 
								fmt2.format(rmsTemp) + " " + fmt2.format(rmsConc) + " " +
								fmt2.format(rmsNo1) + " ");
						double[] bestLoopToNat = bestRMSinLoopGroup(reference, loopSubDir, takeToModel);
						int bestModelToNat = ((int) bestLoopToNat[0]);
						bw.write(fmt2.format(bestLoopToNat[1]) + " " + bestModelToNat + "\n");
						System.out.println(fmt2.format(bestLoopToNat[1]) + " " + bestModelToNat);
					}
				}
				catch (Exception e) {
					// Do nothing
				}
			}
		}
		bw.close();
	}

	
	private static void	printResults(int[] resCount, Protein template) {
		for (int c=0 ; c<resCount.length ; c++) {
			if ((template.residue(c)!=null) && (!template.residue(c).dummy())) {
				System.out.println(c + " " + resCount[c]);				
			}			
		}		
	}
	
	private static double calcRMSonHonigBackbone(Protein loopProt, Protein templateProt) {
		double totRms = 0.0;
		int ntot = 0;
		for (int c=loopProt.firstResidue(); c<=loopProt.lastResidue() ; c++) {
			for (int d=0; d<loopProt.residue(c).atoms().size() ; d++) 
				if (!loopProt.residue(c).atoms().atomAt(d).isHydrogen &&
						loopProt.residue(c).atoms().atomAt(d).isBackbone &&
						!loopProt.residue(c).atoms().atomAt(d).name().equals("CB")) {
					Atom atom = loopProt.residue(c).atoms().atomAt(d);
					Atom atomr = templateProt.residue(c).atoms().findAtomInList(loopProt.residue(c).atoms().atomAt(d).name(),c);
					if (atomr!=null) {
						totRms += (atom.x() - atomr.x())*(atom.x() - atomr.x()) + 
						(atom.y() - atomr.y())*(atom.y() - atomr.y()) +
						(atom.z() - atomr.z())*(atom.z() - atomr.z());
						ntot++;
					}
				}
		}
		return Math.sqrt(totRms/ntot);
	}		

	private static void updateCount(int[] resCount,Protein  loop,Protein  template,double TH,int[] weights) {
		int counter = 0;
		for (int c=loop.firstResidue(); c<=loop.lastResidue() ; c++) {
			Atom atomLoopCA = loop.residue(c).atoms().findAtomInList("CA",c);
			Atom atomLoopCB = loop.residue(c).atoms().findAtomInList("CB",c);
			Atom atomTempCA = template.residue(c).atoms().findAtomInList("CA",c);
			Atom atomTempCB = template.residue(c).atoms().findAtomInList("CB",c);
			if (atomLoopCB!=null) {
				if ((atomLoopCA.distanceFrom(atomTempCA) < TH) && (atomLoopCB.distanceFrom(atomTempCB) < TH)) {
					resCount[c] += weights[counter];
				}
			}
			else {
				if (atomLoopCA.distanceFrom(atomTempCA) < TH) {
					resCount[c] += weights[counter];
				}				
			}
			counter++;
		}		
	}
	
	private static int howManyClosures(String[] logFile) {
		int counter999 = 0;
		for (int c=0 ; c<logFile.length ; c++) {
			if (logFile[c].startsWith("999111111"))
				counter999++;
		}
		return counter999++;
	}

	/**
	 * Returning the {bestModel , bestRMS} in a double array. 
	 * @param template
	 * @param loopSubDir
	 * @param modelsToConsider
	 * @return
	 */
	private static double[] bestRMSinLoopGroup(Protein template, String loopSubDir, int modelsToConsider) {
		int minimalModelNumber = -1;
		String[] logFile = File2StringArray.f2a(loopSubDir+"/log.txt");
		double bestRMS = 99999;
		int bestModel = -1;
		double[] returnAR = {-1,-1};
		if (howManyClosures(logFile)>minimalModelNumber) {
			for (int modelC=0 ; modelC<=modelsToConsider ; modelC++) {
				if ((new File(loopSubDir+"/"+modelC+".pdb.min")).exists()) {
					Protein loop = new Protein(loopSubDir+"/"+modelC+".pdb.min", new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
					double rms = calcRMSonHonigBackbone(loop, template);
					if (rms<bestRMS) {
						bestRMS = rms;
						bestModel = modelC;
					}
				}					
			}
			returnAR[0] = bestModel;
			returnAR[1] = bestRMS; 			 
		}
		else {
			//System.out.println("Best model is: " + loopSubDir+"/"+bestModel+".pdb.min     RMS: " + bestRMS);
		}
		return returnAR;		
	}
	
	private static Protein extractLoopFromTemplate(Protein template , int start , int end) {
		AtomList al = new AtomList();
		for (int c = 0 ; c<template.atoms().size() ; c++) {
			if ((template.atoms().atomAt(c).residueNumber() >= start) &&
					(template.atoms().atomAt(c).residueNumber() <= end)) {
				al.add(new Atom(template.atoms().atomAt(c)));
			}
		}
		return new Protein(al , new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
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
		return (new Protein(tmpprot.atoms().filter(new AtomList.IsDefrostedFilter()).duplicate() , new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS)));
	}	
	
	private static Protein mergeTemplateAndLoop(Protein shortTemplate, Protein loop) {
		AtomList newTemplate = shortTemplate.atoms().duplicate();
		for (int c=0; c<loop.atoms().size() ; c++) {
			Atom atom = newTemplate.findAtomInList(loop.atoms().atomAt(c).name(), 
					loop.atoms().atomAt(c).residueNumber());
			atom.setXYZ(loop.atoms().atomAt(c).x(),
					loop.atoms().atomAt(c).y(),
					loop.atoms().atomAt(c).z());
		}
		Protein newProtein = new Protein(newTemplate, new ResidueExtendedAtoms(ADD_ATOMS));
		RotamericTools.jumble(newProtein);
		for (int c=0; c<newProtein.atoms().size() ; c++) {
			newProtein.atoms().atomAt(c).setChain("A");
		}
		return newProtein;
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
				"Usage java -Xmx600m CalculateLoopConsistencyOnProtein <directory of Loops> \n"+
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));

		loopDirString = getOrderedArgument(args);
		if (loopDirString == null) throw new RuntimeException(errorMessage);
		System.out.println("# Loop dir is: "+loopDirString);

		loopStartString = getOrderedArgument(args);

		loopEndString = getOrderedArgument(args);

		initRandom(999);
	}		
	
}
