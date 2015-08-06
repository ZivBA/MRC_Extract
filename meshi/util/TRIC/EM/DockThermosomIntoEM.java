package meshi.util.TRIC.EM;

import java.io.IOException;

import meshi.applications.TriC.TricAlignment;
import meshi.applications.prediction.GDTcalculator;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.parameters.Residues;
import meshi.util.TRIC.PutUnitInAnyTopPosition;
import meshi.util.dssp.DSSP;
import meshi.util.file.MeshiWriter;

public class DockThermosomIntoEM implements Residues {

	private boolean[] toTake = null;
	
	public DockThermosomIntoEM() {}

	public void setToTake(int[][] keepParts , AtomList fullUnit, DSSP dsspUnit, String validSS) {
		toTake = new boolean[1000];
		char chainID  = fullUnit.atomAt(0).chain().trim().charAt(0);
		boolean[] auxToTake = new boolean[1000];
		for (int c=0 ; c<auxToTake.length ; c++) {
			toTake[c] = false;
			auxToTake[c] = false;
		}
		for (int rangeC=0 ; rangeC<keepParts.length ; rangeC++) {
			for (int c=keepParts[rangeC][0] ; c<=keepParts[rangeC][1] ; c++) {
				auxToTake[c] = true;
			}
		}
		for (int c=0 ; c<toTake.length ; c++) {
			if ((fullUnit.findAtomInList("CA", c)!=null) &&
					auxToTake[c] &&
					((validSS.length()==0) || 
					 ((validSS.indexOf(dsspUnit.SSofRes(c, chainID)) != -1)  /*&& dsspUnit.notInEdge(c, chainID)*/))) {
// SetupF		if ((fullUnit.findAtomInList("CA", c)!=null) &&
//					auxToTake[c] &&
//					((validSS.length()==0) || (validSS.indexOf(dsspUnit.SSofRes(c, 'A')) != -1)) &&
//					dsspUnit.notInEdge(c, 'A') &&
//					(fullUnit.findAtomInList("CA", c).z()>-72.0)) {
				toTake[c] = true;
				System.out.print(c + " ");
			}
		}
		System.out.println();
		TricAlignment tricAlignment = new TricAlignment();
		for (int c=0 ; c<toTake.length ; c++) {
			if (toTake[c]) {
				System.out.print(tricAlignment.getNewResNum('I', c, 'A') + " ");
			}
		}
		System.out.println();
		for (int c=0 ; c<toTake.length ; c++) {
			if (toTake[c]) {
				System.out.print(tricAlignment.getNewResNum('I', c, 'B') + " ");
			}
		}
		System.out.println();
		for (int c=0 ; c<toTake.length ; c++) {
			if (toTake[c]) {
				System.out.print(tricAlignment.getNewResNum('I', c, 'G') + " ");
			}
		}
		System.out.println();
		for (int c=0 ; c<toTake.length ; c++) {
			if (toTake[c]) {
				System.out.print(tricAlignment.getNewResNum('I', c, 'D') + " ");
			}
		}
		System.out.println();
		for (int c=0 ; c<toTake.length ; c++) {
			if (toTake[c]) {
				System.out.print(tricAlignment.getNewResNum('I', c, 'E') + " ");
			}
		}
		System.out.println();
		for (int c=0 ; c<toTake.length ; c++) {
			if (toTake[c]) {
				System.out.print(tricAlignment.getNewResNum('I', c, 'H') + " ");
			}
		}
		System.out.println();
		for (int c=0 ; c<toTake.length ; c++) {
			if (toTake[c]) {
				System.out.print(tricAlignment.getNewResNum('I', c, 'Q') + " ");
			}
		}
		System.out.println();
		for (int c=0 ; c<toTake.length ; c++) {
			if (toTake[c]) {
				System.out.print(tricAlignment.getNewResNum('I', c, 'Z') + " ");
			}
		}
	}
	

	/**
	 * Returns the part of the full unit that are 'true' in the toTake array
	 * @param fullUnit
	 * @param unitID
	 * @return
	 */
	public AtomList takePartOfUnit(AtomList fullUnit, char unitID) {
		TricAlignment tricAlignment = new TricAlignment();
		AtomList sameHomoUnit = new AtomList();
		for (int c=0 ; c<fullUnit.size() ; c++) {
			int resNum = fullUnit.atomAt(c).residueNumber();
			if (toTake[tricAlignment.getNewResNum(unitID, resNum, 'I')]) {
				sameHomoUnit.add(fullUnit.atomAt(c));
			}
		}
		return sameHomoUnit;		
	}	
	
	
	/**
	 * Returns the part of the full unit that is in 'keepParts'
	 * @param fullUnit
	 * @param keepParts - should be according to the Thermosome A chain.
	 * @param unitID
	 * @return
	 */
	public AtomList takePartOfUnit(AtomList fullUnit, int[][] keepParts , char unitID) {
		TricAlignment tricAlignment = new TricAlignment();
		int[][] newKeepParts = new int[keepParts.length][2];
		for (int rangeC=0 ; rangeC<keepParts.length ; rangeC++) {
			newKeepParts[rangeC][0] = tricAlignment.getNewResNum('I', keepParts[rangeC][0], unitID);
			newKeepParts[rangeC][1] = tricAlignment.getNewResNum('I', keepParts[rangeC][1], unitID);
		}
		AtomList sameHomoUnit = new AtomList();
		for (int c=0 ; c<fullUnit.size() ; c++) {
			int resNum = fullUnit.atomAt(c).residueNumber();
			boolean keep = false;
			for (int rangeC=0 ; rangeC<newKeepParts.length ; rangeC++) {
				if ((resNum>=newKeepParts[rangeC][0]) && (resNum<=newKeepParts[rangeC][1])) {
					keep = true;
				}
			}
			if (keep) {
				sameHomoUnit.add(fullUnit.atomAt(c));
			}
		}
		return sameHomoUnit;		
	}
	
	
	/**
	 * 
	 * @param yaoUnit 
	 * @param keepParts - should be according to the Thermosome A chain.
	 * @return
	 */
	public AtomList dockSameHomologyUnit(AtomList yaoUnit, int[][] keepParts, double cutoff) {
		String chainID = yaoUnit.atomAt(0).chain();
		AtomList sameHomoUnitFull = new AtomList("Atemp_HM_"+chainID.charAt(0)+".pdb");
		AtomList sameHomoUnit = takePartOfUnit(sameHomoUnitFull, keepParts, chainID.charAt(0));
		GDTcalculator.alignBySubset(yaoUnit, sameHomoUnit, cutoff);
		AtomList returnSameHomoUnit = new AtomList("Atemp_HM_"+chainID.charAt(0)+".pdb");
		GDTcalculator.alignBySubset(sameHomoUnit, returnSameHomoUnit, 0.0001);
		returnSameHomoUnit.setChain(chainID);
		return returnSameHomoUnit;
	}
	
	public void fixResNumbering(AtomList atoms) {
		String chainID = atoms.atomAt(0).chain();
		TricAlignment tricAlignment = new TricAlignment();
		String seq = tricAlignment.getAlignment(chainID).replace("-", "");
		String yaoSeq = atoms.getSequence();
		int index1 = seq.indexOf(yaoSeq);
		if (index1 == -1) {
			throw new RuntimeException("No Match");
		}
		if (index1 > 0) { 
			for (int c=0 ; c<atoms.size() ; c++) {
				if (chainID.equals("E")&&(atoms.atomAt(c).residueNumber()>320)) {
					atoms.atomAt(c).setResidueNumber(atoms.atomAt(c).residueNumber() + index1-5);
				}
				else {
					atoms.atomAt(c).setResidueNumber(atoms.atomAt(c).residueNumber() + index1);
				}
			}				
		}
	}
	
	
	
	public static void main(String[] args) {
		int zvl = ALA;
		String yaoPositions = "QDEHAZBG";
		String units = "ABGDEHQZ";
// Setup1,2		int[][] keepParts = {{17 , 144}, {403 , 520}};  // Always in the thermoA frame.
// Setup3		int[][] keepParts = {{260 , 285}};  // Always in the thermoA frame.
// Setup4		int[][] keepParts = {{147 , 215}, {366 , 404}};  // Always in the thermoA frame.
// SetupAB		int[][] keepParts = {{215 , 244} , {274 , 366}};  // Always in the thermoA frame.
// SetupCD		int[][] keepParts = {{147 , 183} , {196 , 215} , {366 , 403}};  // Always in the thermoA frame.
// SetupEF		int[][] keepParts = {{17 , 144}, {403 , 520}};  // Always in the thermoA frame.
// SetupGH		int[][] keepParts = {{252 , 257}, {260 , 275}};  // Always in the thermoA frame.
		int[][] keepParts = {{17 , 144}, {403 , 520}};  // Always in the thermoA frame.
//		int[][] keepPartsFull = {{1 , 525}};  // Always in the thermoA frame.
		String fileNamePrefix = "SetupE";
		String fileNamePrefixWrite = "SetupE";
//		AtomList yaoFull = new AtomList(args[0].trim());
//		System.out.println("Initial scaffold: " + args[0].trim());
		DockThermosomIntoEM dock = new DockThermosomIntoEM();
		
		
		// Getting the fitSet
		DSSP thermoAdssp = new DSSP("1A6D_chainA.dssp");
		AtomList thermoA = new AtomList("1A6D_chainA.pdb");
// SetupF		AtomList thermoA = new AtomList("Setup2_pos_0_full.pdb");
		dock.setToTake(keepParts, thermoA, thermoAdssp, "HE");
/* Setup1,2		
  		dock.toTake[37]=dock.toTake[38]=dock.toTake[487]=dock.toTake[488]=false;
		dock.toTake[410]=dock.toTake[441]=dock.toTake[442]=dock.toTake[443]=true;
Setup1,2 */
/* Setup4				dock.toTake[158]=dock.toTake[380]=dock.toTake[381]=true;
		for (int c=0 ; c<dock.toTake.length ; c++) {
			dock.toTake[c] = false;
		}
		for (int c=184 ; c<192 ; c++) {
			dock.toTake[c] = true;
		}
Setup4	*/		
/* Setup5
		for (int c=0 ; c<dock.toTake.length ; c++) {
			dock.toTake[c] = false;
		}
		dock.toTake[295]=dock.toTake[296]=true;
Setup5	*/
/* SetupAB 
//		dock.toTake[329]=dock.toTake[330]=dock.toTake[350]=dock.toTake[351]=dock.toTake[352]=dock.toTake[355]=dock.toTake[356]=false;
//		dock.toTake[286]=dock.toTake[291]=dock.toTake[311]=true;
		for (int c=0 ; c<dock.toTake.length ; c++) {
			dock.toTake[c]=false;
		}
		for (int c=296 ; c<=304 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=315 ; c<=324 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=287 ; c<=290 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=308 ; c<=311 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=230 ; c<=239 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=341 ; c<=349 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=357 ; c<=364 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=283 ; c<=286 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=305 ; c<=307 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=312 ; c<=314 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=294 ; c<=295 ; c++) {
			dock.toTake[c]=true;
		}
SetupAB */
/* SetupCD 	
//		dock.toTake[197]=dock.toTake[214]=dock.toTake[215]=dock.toTake[401]=dock.toTake[402]=dock.toTake[403]=false;
//		dock.toTake[202]=dock.toTake[210]=dock.toTake[374]=dock.toTake[380]=true;
		for (int c=0 ; c<dock.toTake.length ; c++) {
			dock.toTake[c]=false;
		}
		for (int c=147 ; c<=159 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=167 ; c<=182 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=380 ; c<=402 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=197 ; c<=202 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=210 ; c<=213 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=368 ; c<=374 ; c++) {
			dock.toTake[c]=true;
		}		
SetupCD */
/* SetupEF */ 
//		dock.toTake[65]=dock.toTake[66]=dock.toTake[67]=dock.toTake[68]=dock.toTake[69]=false;
//		dock.toTake[77]=dock.toTake[78]=dock.toTake[79]=dock.toTake[80]=dock.toTake[81]=false;
//  	dock.toTake[37]=dock.toTake[38]=dock.toTake[487]=dock.toTake[488]=false;
//		dock.toTake[410]=dock.toTake[441]=dock.toTake[442]=dock.toTake[443]=true;
		for (int c=0 ; c<dock.toTake.length ; c++) {
			dock.toTake[c]=false;
		}
		for (int c=21 ; c<=37 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=95 ; c<=116 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=120 ; c<=140 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=410 ; c<=424 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=429 ; c<=452 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=456 ; c<=468 ; c++) {
			dock.toTake[c]=true;
		}
		for (int c=496 ; c<=513 ; c++) {
			dock.toTake[c]=true;
		}
/*SetupEF */
		System.out.println();
		int trueCounter=0;
		for (int c=0 ; c<dock.toTake.length ; c++) {
			if (dock.toTake[c]) {
				System.out.print(c + " ");
				trueCounter++;
			}
			if (trueCounter%25==0) {
				System.out.println();
			}
		}
		System.out.println("\n" + trueCounter);
		thermoA.setChain("I");
		
		
//		// Getting thermoA in the 8 positions.
//		for (int yaoC = 0 ; yaoC<8 ; yaoC++) {
//			Atom.resetNumberOfAtoms();
//			String yaoID = yaoPositions.substring(yaoC,yaoC+1);
//			dock.fixResNumbering(yaoFull.chainFilter(yaoID));
//			AtomList sameHomoUnit = dock.dockSameHomologyUnit(yaoFull.chainFilter(yaoID), keepParts, 0.75);
//			PutUnitInAnyTopPosition.alignByEquatorialDomains(sameHomoUnit, yaoID.charAt(0), thermoA, 'I');
//			try {
//				thermoA.print(new MeshiWriter(fileNamePrefixWrite+"_pos_"+yaoC+"_full.pdb"));
//			} catch (IOException e) {
//				throw new RuntimeException("Could not write file.");
//			}				
//		}
		
		
		// Putting the different units in different positions
		// Looping on positions.
		for (int yaoC = 0 ; yaoC<8 ; yaoC++) {
			Atom.resetNumberOfAtoms();
			AtomList thermoAinPos = new AtomList(fileNamePrefix+"_pos_"+yaoC+"_full.minMap.pdb");
			// Looping on units
			for (int unitC = 0 ; unitC<8 ; unitC++) {
				Atom.resetNumberOfAtoms();
				AtomList anyUnitFull = new AtomList("Atemp_HM_"+units.charAt(unitC)+".scwrl.pdb");
				PutUnitInAnyTopPosition.alignByEquatorialDomains(thermoAinPos, 'I', anyUnitFull, units.charAt(unitC));
				AtomList anyUnit = dock.takePartOfUnit(anyUnitFull, units.charAt(unitC));
				anyUnit.setChain(units.charAt(unitC)+"");
				try {
					anyUnit.print(new MeshiWriter(fileNamePrefixWrite+"_"+yaoC+"_"+units.charAt(unitC)+".pdb"));
				} catch (IOException e) {
					throw new RuntimeException("Could not write file.");
				}				
			}
			Atom.resetNumberOfAtoms();
			AtomList anyUnitFull = new AtomList("1A6D_chainA.scwrl.pdb");
			PutUnitInAnyTopPosition.alignByEquatorialDomains(thermoAinPos, 'I', anyUnitFull, 'I');
			AtomList anyUnit = dock.takePartOfUnit(anyUnitFull, 'I');
			anyUnit.setChain("I");
			try {
				anyUnit.print(new MeshiWriter(fileNamePrefixWrite+"_"+yaoC+"_I.pdb"));
			} catch (IOException e) {
				throw new RuntimeException("Could not write file.");
			}							
		}
		
		
		
		
	}
	
	
}




//// Getting thermoA in the 8 positions.
//AtomList newYao = new AtomList();
//TricAlignment tricAlignment = new TricAlignment();
//for (int yaoC = 0 ; yaoC<8 ; yaoC++) {
//	Atom.resetNumberOfAtoms();
//	String yaoID = yaoPositions.substring(yaoC,yaoC+1);
//	dock.fixResNumbering(yaoFull.chainFilter(yaoID));
//	for (int c=0 ; c<dock.toTake.length ; c++) {
//		if (dock.toTake[c]) {
//			Atom atom  = yaoFull.findAtomInListReturningAtom("CA", yaoID, tricAlignment.getNewResNum('I', c, yaoID.charAt(0)));
//			atom.addToZ(0);
//			newYao.add(atom);
//		}
//	}
//}
//try {
//	newYao.print(new MeshiWriter(args[0].trim()+".justH.pdb"));
//} catch (IOException e) {
//	throw new RuntimeException("Could not write file.");
//}				
//System.exit(0);