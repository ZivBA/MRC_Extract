package meshi.applications.loopBuilding.applications;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Vector;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.molecularElements.Atom;
import meshi.molecularElements.Protein;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.MinimizerException;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;

public class FindGapsInModel extends MeshiProgram implements Residues, AtomTypes {

	private static String alignmentFileName = null;  
	
	private static String modelFileName = null;

	private static String outFileName = null;

	private static int residuesAroundGap = -1;

	public static void main(String[] args) throws MinimizerException, LineSearchException{
		init(args); 
		
		Protein model = null; 
		int firstResInQuery = -1;
		String queryAlignment = "";
		int resNumCounterQ=0;
		Vector<Integer> badRes = new Vector<Integer>();


		
		// Reading the alignment file
		model = new ExtendedAtomsProtein(modelFileName,DO_NOT_ADD_ATOMS);
		String[] alignmentStrings = File2StringArray.f2a(alignmentFileName);
		firstResInQuery = (new Integer(alignmentStrings[1])).intValue();
		System.out.println("Alignment in query starts at: " + firstResInQuery);
		queryAlignment = alignmentStrings[3].trim();
		System.out.println("QUERY " + queryAlignment);

		
		// Going over the sequence to find gaps
		resNumCounterQ = firstResInQuery-1;
		Atom lastViewedCA = null;
		for (int position=0; position<queryAlignment.length() ; position++) {
			if (queryAlignment.charAt(position)!='-') {
				resNumCounterQ++;
				if ((model.residue(resNumCounterQ)!= null) && (model.residue(resNumCounterQ).ca()!= null)) {
					if (lastViewedCA!=null) {
						if (Math.abs(model.residue(resNumCounterQ).ca().distanceFrom(lastViewedCA) - 3.75)>0.25) // A gap here
							badRes.add(new Integer(resNumCounterQ));
					}
					lastViewedCA = model.residue(resNumCounterQ).ca();
				}
				else { // Missing residue
					badRes.add(new Integer(resNumCounterQ));
					lastViewedCA = null;
				}
			}
		}
		
		
		
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFileName));
			int gapStart = -999;
			int gapEnd = -999;
			for (int cc=0 ; cc<badRes.size() ; cc++) {
				if (badRes.get(cc).intValue()!=(gapEnd+1)) {
					if (gapStart!=-999) {
						if (gapStart!=firstResInQuery) {
							bw.write((gapStart-residuesAroundGap) + " " + (gapEnd+residuesAroundGap) + " 999 \n");
							System.out.println("BAD: " + (gapStart-residuesAroundGap) + " " + (gapEnd+residuesAroundGap) + " 999 ");
						}
						else {
							bw.write(gapStart + " " + (gapEnd+residuesAroundGap) + " 999 \n");
							System.out.println("BAD: " + gapStart + " " + (gapEnd+residuesAroundGap) + " 999 ");							
						}
					}
					gapStart = badRes.get(cc).intValue();
					gapEnd = gapStart;
				}
				else {
					gapEnd++;
				}
			}
			if (gapStart!=-999) {
				if (gapStart==firstResInQuery) {
					bw.write(gapStart + " " + (gapEnd+residuesAroundGap) + " 999 \n");
					System.out.println("BAD: " + gapStart + " " + (gapEnd+residuesAroundGap) + " 999 ");							
				}
				else if (gapEnd==resNumCounterQ) {
					System.out.println("BAD: " + (gapStart-residuesAroundGap) + " " + gapEnd + " 999 ");
					bw.write((gapStart-residuesAroundGap) + " " + gapEnd + " 999 \n");					
				}
				else {
					System.out.println("BAD: " + (gapStart-residuesAroundGap) + " " + (gapEnd+residuesAroundGap) + " 999 ");
					bw.write((gapStart-residuesAroundGap) + " " + (gapEnd+residuesAroundGap) + " 999 \n");
				}
			}
			bw.close();			
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
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
				"Usage java -Xmx600m FindGapsInModel <alignment file name> <model filename> <output filename> <take buffer around gap>\n"+
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));

		alignmentFileName = getOrderedArgument(args);
		if (alignmentFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Alignment file name is: "+alignmentFileName);

		modelFileName = getOrderedArgument(args);
		if (modelFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Completeing: "+modelFileName);
		
		outFileName = getOrderedArgument(args);
		if (outFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Output to: "+outFileName);

		String tmp = getOrderedArgument(args);
		if (tmp == null) throw new RuntimeException(errorMessage);
		residuesAroundGap = (new Integer(tmp.trim())).intValue();
		System.out.println("# Residues to take around gap: "+residuesAroundGap);
		
		initRandom(999);
	}	

} // Of AddMissingResidues
