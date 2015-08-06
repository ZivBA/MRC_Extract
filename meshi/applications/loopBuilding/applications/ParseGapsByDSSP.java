package meshi.applications.loopBuilding.applications;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Vector;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.dssp.DSSP;

public class ParseGapsByDSSP extends MeshiProgram implements Residues, AtomTypes {


	private static String dsspFileName = null;
	private static String initialFileName = null;
	private static String gapsFileName = null;
	private static int minimalHelix = 7;
	private static int minimalStrand = 2;
	private static int removeHelix = 2;
	private static int removeStrand = 1;
	private static final int MIN_LOOP_SIZE = 4;   
	private static final int STANDARD_LOOP_SIZE = 8;   
	private static final int TRIM_LOOP_END = 2;   
	

	
	public static void main(String[] args) {
		init(args); 
		
		if ((removeHelix<1) || (removeStrand<1))
			throw new RuntimeException("The removal must be 1 or more.");

		Vector<String>  output_superLoops = new Vector<String>();
		Vector<String>  output_superLoops_detail = new Vector<String>();
		Vector<int[]>  output_superLoops_limits = new Vector<int[]>();
		DSSP dssp = new DSSP(dsspFileName);	
		Protein model = new ExtendedAtomsProtein(initialFileName,DO_NOT_ADD_ATOMS);
		int modelFirstRes = model.firstResidue();
		int modelLastRes = ((Residue) model.residues().last()).number;
		
		// Finding missing residues
		int firstRes = modelFirstRes;
		for (int gapRes=modelFirstRes; gapRes<=modelLastRes+1 ; gapRes++) {
			if (!allKosher(model,gapRes,gapRes)) {
				System.out.println("chain [" + firstRes + "," + (gapRes-1) + "]:");
				int lastRes = (gapRes-1);

				// Finding the secondary structures
				Vector<int[]> sss = new Vector<int[]>();
				char lastSS = 'X';
				for (int res=firstRes; res<=lastRes ; res++) {
					if ((dssp.SSofRes(res,' ')=='E')||(dssp.SSofRes(res,' ')=='H')) {
						if (dssp.SSofRes(res,' ')!=lastSS)  {
							int startSS = res;
							int endSS = res;
							for (endSS = res; (dssp.SSofRes(res,' ')==dssp.SSofRes(endSS,' ')) && (endSS<=lastRes) ; endSS++) {}
							endSS--;
							if (((dssp.SSofRes(res,' ')=='E') && ((endSS-startSS+1) > minimalStrand)) || 
									((dssp.SSofRes(res,' ')=='H') && ((endSS-startSS+1) > minimalHelix))) {
								System.out.println("Added SS of " + dssp.SSofRes(res,' ') + " at " + "[" + startSS +
										"," + endSS + "].");
								if ((dssp.SSofRes(res,' ')=='E')) {
									int[] tmp = {startSS+removeStrand,endSS-removeStrand};
									sss.add(tmp);
								}
								if ((dssp.SSofRes(res,' ')=='H')) {
									int[] tmp = {startSS+removeHelix,endSS-removeHelix};
									sss.add(tmp);
								}
							}
							else {
								System.out.println("The SS of " + dssp.SSofRes(res,' ') + " at " + "[" + startSS +
										"," + endSS + "] is too small - ignored.");
							}					
							res = endSS;
						}
					}
					lastSS = dssp.SSofRes(res,' ');
				}
				
				// Finding the GAPS
				// Is there an opening gap?
				if ((sss.firstElement()[0]-firstRes)<(MIN_LOOP_SIZE/2)) {
					System.out.println("Unable to refine ["+firstRes+","+(sss.firstElement()[0]-1)+"] - too short");			
				} else {
					System.out.println("GAP ["+firstRes+","+(sss.firstElement()[0]-1)+"] - free end");
					doOutputStrings(firstRes, (sss.firstElement()[0]-1), true,
							output_superLoops, output_superLoops_detail, output_superLoops_limits);					
				}
				// Doing all other loops
				for (int ss=0; ss<sss.size()-1 ; ss++) {
					if ((sss.get(ss+1)[0]-sss.get(ss)[1]-1)<MIN_LOOP_SIZE) {
						System.out.println("Unable to refine ["+(sss.get(ss)[1]+1)+","+(sss.get(ss+1)[0]-1)+"] - too short");			
					} else {
						System.out.println("GAP ["+(sss.get(ss)[1]+1)+","+(sss.get(ss+1)[0]-1)+"]");
						doOutputStrings((sss.get(ss)[1]+1),(sss.get(ss+1)[0]-1), false,
								output_superLoops, output_superLoops_detail, output_superLoops_limits);					
					}					
				}
				// Is there an ending gap?
				if ((lastRes - sss.lastElement()[1])<(MIN_LOOP_SIZE/2)) {
					System.out.println("Unable to refine ["+(sss.lastElement()[1]+1)+","+lastRes+"] - too short");			
				} else {
					System.out.println("GAP ["+(sss.lastElement()[1]+1)+","+lastRes+"] - free end");
					doOutputStrings((sss.lastElement()[1]+1),lastRes, true,
							output_superLoops, output_superLoops_detail, output_superLoops_limits);					
				}
				
				
				for ( ; (gapRes<=lastRes+1) && !allKosher(model,gapRes,gapRes); gapRes++) {}
				firstRes = gapRes;
				gapRes--;
			}
		}
		
		// Printing the output
		for (int c=0 ; c<output_superLoops.size() ; c++) {
			System.out.println("Title: " + output_superLoops_limits.get(c)[0] + " " + output_superLoops_limits.get(c)[1]);
			System.out.print("Master line: " + output_superLoops.get(c));
			System.out.println("Details:");
			System.out.print(output_superLoops_detail.get(c));
		}
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(gapsFileName));
			for (int c=0 ; c<output_superLoops.size(); c++) {
				bw.write(output_superLoops.get(c));
			}
			bw.close();			
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}			
		try{
			for (int c=0 ; c<output_superLoops.size(); c++) {
				BufferedWriter bw = new BufferedWriter(new FileWriter(gapsFileName+"."+output_superLoops_limits.get(c)[0]+"_"+
							output_superLoops_limits.get(c)[1]+".txt"));
				bw.write(output_superLoops_detail.get(c));
				bw.close();			
			}
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}			
	}
	
	private static boolean allKosher(Protein prot, int start, int end) {
		for (int c=start ; c<=end ; c++)
			if ((prot.residue(c)==null) || (prot.residue(c).ca()==null))
				return false;
		return true;
	}
	
	private static void doOutputStrings(int resStart, int resEnd, boolean freeEnd,
			Vector<String>  output_superLoops, Vector<String>  output_superLoops_detail,
			Vector<int[]>  output_superLoops_limits) {
		int[] tmp = {resStart,resEnd};
		output_superLoops_limits.add(tmp);
		String out ="";
		if (freeEnd) {
			out += (resStart + " " + resEnd + " 0 " + resStart + " " + resEnd + " 999\n");			
			output_superLoops.add(resStart + " " + resEnd + " 1 999\n");
		}
		else {
			if ((resEnd-resStart+1)<=STANDARD_LOOP_SIZE) {
				out += (resStart + " " + resEnd + " 0 " + resStart + " " + resEnd + " 999\n");
				output_superLoops.add(resStart + " " + resEnd + " 1 999\n");
			}
			else {
				output_superLoops.add(resStart + " " + resEnd + " 0 999\n");
				int manner = -1;
				int start = resStart;
				int end = resEnd;
				while (end>start) {
					if (manner==-1) {
						int END = Math.min(start+STANDARD_LOOP_SIZE-1,end+TRIM_LOOP_END);
						out += (start + " " + END + " -1 " + start + " " + (END-TRIM_LOOP_END) + " 999\n");
						start = (END-TRIM_LOOP_END+1);
					}
					else {
						int START = Math.max(end-STANDARD_LOOP_SIZE+1, start-TRIM_LOOP_END);
						out += (START + " " + end + " 1 " + (START+TRIM_LOOP_END) + " " + end + " 999\n");
						end = (START+TRIM_LOOP_END-1);
					}
					manner = -manner;
				}
			}
		}
		output_superLoops_detail.add(out);
	}


	/** ================================= init =========================================
	 *
	 *A static function for parsing of the command line arguments and assigning the 
	 *variables commandsFileName, modelFileName and randomNumberSeed with the right inputs. Note that this
	 *static method is using parsing functions such as getOrderedArguments that are defined in MeshiProgram
	 *that MinimizeProtein inherits.
	 **/

	protected static void init(String[] args) {

		String errorMessage = ("\n                  ******************\n"+
				"Usage java -Xmx300m MinimizeBatchOfCoarseLoops <commands file name> <file with list of pdbs> <file with list of phipsi data> <ref pdb file name> " +
				"<loop starting resisue> <loop ending residue> <Wev> <Whb> <Wtorval> <Wtether> <output PDB extension string>\n"+
		"                    ******************\n");

		int zvl = ALA; // force the reading of "meshi.parameters.Residues"
		zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"

		initRandom(999);
		
		dsspFileName = getOrderedArgument(args);
		if (dsspFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# The DSSP file name is: "+dsspFileName);

		initialFileName = getOrderedArgument(args);
		if (initialFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# The initial file name model: "+initialFileName);
		
		gapsFileName = getOrderedArgument(args);
		if (gapsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# The gaps file name: "+gapsFileName);

		String tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		minimalHelix = (new Integer(tmpString)).intValue();
		System.out.println("# Minimial HELIX length: " + minimalHelix);
		
		tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		minimalStrand = (new Integer(tmpString)).intValue();
		System.out.println("# Minimial STRAND length: " + minimalStrand);
		
		tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		removeHelix = (new Integer(tmpString)).intValue();
		System.out.println("# Trim HELIX by: " + removeHelix);
		
		tmpString = getOrderedArgument(args);
		if (tmpString== null) throw new RuntimeException(errorMessage);
		removeStrand = (new Integer(tmpString)).intValue();
		System.out.println("# Minimial STRAND length: " + removeStrand);
		
	}	

} // Of AddMissingResidues
