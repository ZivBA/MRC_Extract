package meshi.util.external;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.StringTokenizer;

import meshi.applications.TriC.TricAlignment;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionList;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;
import meshi.util.rotamericTools.RotamericTools;

public class TINKERandMESHIcomparison implements Residues {
	
	private final double posRestrainSpringConst = 100.0;	
	private final double posRestrainTolerance = 0.5;   
	private final double torRestrainSpringConst = 1.0;
	String inactiveString = "INACTIVE ";
	String restrainString = "RESTRAIN-POSITION ";
	String restrainTorString = "RESTRAIN-TORSION ";
	/**
	 * How to convert the indices of the atoms in the MESHI AtomList to the indices in the XYZ file. 
	 * (-1) if no correspondence exists.
	 */
	private int[] meshi2tinker = null;
	/**
	 * How to convert the indices of the atoms in the XYZ file to the indices in the MESHI AtomList. 
	 * (-1) if no correspondence exists.
	 */
	private int[] tinker2meshi = null; 
	private AtomList alMeshi;  // Atom list of the PDB file
	private String tinkerXYZ = null; // The tinker's XYZ file
	/**
	 * A record of the XYZ file atoms that are inactive.
	 */
	private boolean[] inactive = null;
	/**
	 * connected[c] gives an array of the connected atoms to XYZ atom 'c'. All numbers in XYZ file indices.
	 */
	private int[][] connected = null;

	
	public TINKERandMESHIcomparison(String tinkerXYZ, String meshiPDB) {
		this.tinkerXYZ = tinkerXYZ; 
		alMeshi = new AtomList(meshiPDB);
		meshi2tinker = new int[alMeshi.size()];
		for (int c=0; c<alMeshi.size() ; c++) {
			meshi2tinker[c] = -1;
		}
		String[] lines = File2StringArray.f2a(tinkerXYZ);
		int numberOfTinkerAtoms = Integer.parseInt(lines[0].trim());
		tinker2meshi = new int[numberOfTinkerAtoms+1];
		connected = new int[numberOfTinkerAtoms+1][];
		for (int c=0; c<tinker2meshi.length ; c++) {
			tinker2meshi[c] = -1;
		}
		for (int c=1; c<=numberOfTinkerAtoms ; c++) {
			StringTokenizer st = new StringTokenizer(lines[c]);
			if (Integer.parseInt(st.nextToken())!=c){
				throw new RuntimeException("Atom numbers do not match.");
			}
			st.nextToken(); //The atom type 
			double x = Double.parseDouble(st.nextToken());
			double y = Double.parseDouble(st.nextToken());
			double z = Double.parseDouble(st.nextToken());
			double minDis = Double.MAX_VALUE;
			int minInd = -1;
			for (int ac=0; ac<alMeshi.size() ; ac++) {
				double ax = alMeshi.atomAt(ac).x();
				double ay = alMeshi.atomAt(ac).y();
				double az = alMeshi.atomAt(ac).z();
				double dis = (ax-x)*(ax-x) + (ay-y)*(ay-y) + (az-z)*(az-z);
				if (dis<minDis) {
					minDis = dis;
					minInd = ac;
				}				
			}
			if (minDis<0.03) {
//				System.out.println("Match:\n" + lines[c] + "\n" + alMeshi.atomAt(minInd));
				tinker2meshi[c] = minInd;
				meshi2tinker[minInd] = c;
			}
			st.nextToken(); //The atom OPLS number 
			connected[c] = new int[st.countTokens()];
			for (int cc=0 ; st.hasMoreTokens() ; cc++) {
				connected[c][cc] = Integer.parseInt(st.nextToken());
			}			
		}
		// Setting all the atoms to active state
		inactive = new boolean[tinker2meshi.length];
		for (int c=0 ; c<tinker2meshi.length ; c++) {
			inactive[c] = false;
		}
	}
	
	
	public int tinker2meshi(int tinkerInd) {
		return tinker2meshi[tinkerInd];
	}
	
	public int meshi2tinker(int meshiInd) {
		return meshi2tinker[meshiInd];
	}
	
	/*
	 * This method will inactive the thermosom superposition atoms
	 */
	public void setInactiveSuperPositionAtoms(char unitName) {
		TricAlignment tric = new TricAlignment();		
		alMeshi.defrost();
		inactive = new boolean[tinker2meshi.length];
		for (int c=0 ; c<tinker2meshi.length ; c++) {
			inactive[c] = false;
		}
    	int[][] ThermoA = {{24,34},
    			{97,113},
    			{123,139}}; 
    	for (int res=24 ; res<=34 ; res++) {
    		int resInUnit = tric.getNewResNum('I', res, unitName);
    		int indInList = alMeshi.findAtomInList("CA", "A", resInUnit);
    		alMeshi.atomAt(indInList).freeze();
    		inactive[meshi2tinker[indInList]] = true;
    	}
    	for (int res=97 ; res<=113 ; res++) {
    		int resInUnit = tric.getNewResNum('I', res, unitName);
    		int indInList = alMeshi.findAtomInList("CA", "A", resInUnit);
    		alMeshi.atomAt(indInList).freeze();
    		inactive[meshi2tinker[indInList]] = true;
    	}
    	for (int res=123 ; res<=139 ; res++) {
    		int resInUnit = tric.getNewResNum('I', res, unitName);
    		int indInList = alMeshi.findAtomInList("CA", "A", resInUnit);
    		alMeshi.atomAt(indInList).freeze();
    		inactive[meshi2tinker[indInList]] = true;
    	}
	}
	
	
	/*
	 * This method will also freeze the inactive atoms in 'alMeshi'
	 */
	public void setInactiveFarFromInterface(double disCutOff) {
		alMeshi.defrost();
		inactive = new boolean[tinker2meshi.length];
		for (int c=0 ; c<tinker2meshi.length ; c++) {
			inactive[c] = false;
		}
		for (int ac=0; ac<alMeshi.size() ; ac++) {
			double disFromInterface = Double.MAX_VALUE;
			for (int ac1=0; ac1<alMeshi.size() ; ac1++) {
				double x = alMeshi.atomAt(ac).x();
				double y = alMeshi.atomAt(ac).y();
				double z = alMeshi.atomAt(ac).z();
				double ax = alMeshi.atomAt(ac1).x();
				double ay = alMeshi.atomAt(ac1).y();
				double az = alMeshi.atomAt(ac1).z();
				double dis = Math.sqrt((ax-x)*(ax-x) + (ay-y)*(ay-y) + (az-z)*(az-z));
				if (!alMeshi.atomAt(ac).chain().equals(alMeshi.atomAt(ac1).chain())) {
					if (dis<disFromInterface) {
						disFromInterface = dis;
					}					
				}
			}
			if (disFromInterface>disCutOff) {
				inactive[meshi2tinker(ac)] = true;
				alMeshi.atomAt(ac).freeze();
				for (int cc=0 ; cc<connected[meshi2tinker(ac)].length ; cc++) {
					if (tinker2meshi[connected[meshi2tinker(ac)][cc]]==-1) {
						inactive[connected[meshi2tinker(ac)][cc]] = true;
					}
				}
			}
		}
	}
	

	/*
	 * This method will copy the inactivation data from the other object
	 */
	public void setInactiveFarFromInterface(TINKERandMESHIcomparison othertink) {
		alMeshi.defrost();
		inactive = new boolean[tinker2meshi.length];
		for (int c=0 ; c<tinker2meshi.length ; c++) {
			inactive[c] = othertink.inactive[c];
		}
		for (int ac=0; ac<alMeshi.size() ; ac++) {
			if (inactive[meshi2tinker(ac)]) {
				alMeshi.atomAt(ac).freeze();
			}
		}		
	}

	
	
	public String getInactiveString() {
		String outString = "";
		for (int c=1 ; c<inactive.length ; c++) {
			if (inactive[c]) {
				int lastInactive=c;
				for ( ; ((lastInactive+1)<inactive.length) && inactive[lastInactive+1] ; lastInactive++);
				if (lastInactive==c) 
					outString += (inactiveString + c + "\n"); 
				else 
					outString += (inactiveString + "-" + c + " " + lastInactive + "\n");
				c = lastInactive;
			}
		}
		return outString;		
	}
	
	/**
	 * Returning the restraining stgring for backbone active atoms.
	 */	
	public String getRestrainString() {
		String outString = "";
		for (int c=1 ; c<inactive.length ; c++) {
			if (!inactive[c]) {
				int ac = tinker2meshi(c);
				if ((ac>-1) && alMeshi.atomAt(ac).isBackbone && !alMeshi.atomAt(ac).isHydrogen) {
					outString += (restrainString + c + " " + 
							alMeshi.atomAt(ac).x() + " " + 
							alMeshi.atomAt(ac).y() + " " +
							alMeshi.atomAt(ac).z() + " " +
							posRestrainSpringConst + " " + posRestrainTolerance + "\n");
				}				
			}
		}
		return outString;		
	}
	
	public String getTorsionRestrainString(double[][] pp) {
		String outString = "";
		Protein meshiProtein = ComplexMESHIconversion.complex2meshi(alMeshi);
		DistanceMatrix	distanceMatrix = new DistanceMatrix(meshiProtein.atoms(), 2.0, 1.0, 4);
		TorsionList modelChi1 = (TorsionList) TorsionList.createTorsionList(meshiProtein,distanceMatrix).filter(new TorsionList.FilterChi1() , new TorsionList());
		TorsionList modelChi2 = (TorsionList) TorsionList.createTorsionList(meshiProtein,distanceMatrix).filter(new TorsionList.FilterChi2() , new TorsionList());
		TorsionList modelChi3 = (TorsionList) TorsionList.createTorsionList(meshiProtein,distanceMatrix).filter(new TorsionList.FilterChi3() , new TorsionList());
		TorsionList modelChi4 = (TorsionList) TorsionList.createTorsionList(meshiProtein,distanceMatrix).filter(new TorsionList.FilterChi4() , new TorsionList());
		for (int chiCounter = 0 ; chiCounter<modelChi1.size() ; chiCounter++ ) {
			Torsion tor = modelChi1.torsionAt(chiCounter);
			if (!tor.frozen()) {
				int ind1 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom1)); 
				int ind2 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom2)); 
				int ind3 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom3)); 
				int ind4 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom4)); 
				int torsionNum = tor.getTorsionResNum();
				double weight = torRestrainSpringConst * (10.0*Math.PI/180.0) / pp[torsionNum][8];
				outString += (restrainTorString + ind1 + " " + ind2 + " " + ind3 + " " + ind4 + " " + (weight*Math.PI/180.0) + " " + 
						(pp[torsionNum][4]*180.0/Math.PI - 1e-4) + " " + (pp[torsionNum][4]*180.0/Math.PI + 1e-4) + " " + "\n");
			}
		}
		for (int chiCounter = 0 ; chiCounter<modelChi2.size() ; chiCounter++ ) {
			Torsion tor = modelChi2.torsionAt(chiCounter);
			if (!tor.frozen()) {
				int ind1 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom1)); 
				int ind2 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom2)); 
				int ind3 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom3)); 
				int ind4 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom4)); 
				int torsionNum = tor.getTorsionResNum();
				double weight = torRestrainSpringConst * (10.0*Math.PI/180.0) / pp[torsionNum][9];
				outString += (restrainTorString + ind1 + " " + ind2 + " " + ind3 + " " + ind4 + " " + (weight*Math.PI/180.0) + " " + 
						(pp[torsionNum][5]*180.0/Math.PI - 1e-4) + " " + (pp[torsionNum][5]*180.0/Math.PI + 1e-4) + " " + "\n");
			}
		}
		for (int chiCounter = 0 ; chiCounter<modelChi3.size() ; chiCounter++ ) {
			Torsion tor = modelChi3.torsionAt(chiCounter);
			if (!tor.frozen()) {
				int ind1 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom1)); 
				int ind2 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom2)); 
				int ind3 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom3)); 
				int ind4 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom4)); 
				int torsionNum = tor.getTorsionResNum();
				double weight = torRestrainSpringConst * (10.0*Math.PI/180.0) / pp[torsionNum][10];
				outString += (restrainTorString + ind1 + " " + ind2 + " " + ind3 + " " + ind4 + " " + (weight*Math.PI/180.0) + " " + 
						(pp[torsionNum][6]*180.0/Math.PI - 1e-4) + " " + (pp[torsionNum][6]*180.0/Math.PI + 1e-4) + " " + "\n");
			}
		}
		for (int chiCounter = 0 ; chiCounter<modelChi4.size() ; chiCounter++ ) {
			Torsion tor = modelChi4.torsionAt(chiCounter);
			if (!tor.frozen()) {
				int ind1 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom1)); 
				int ind2 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom2)); 
				int ind3 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom3)); 
				int ind4 = meshi2tinker(ComplexMESHIconversion.findAtomInOrigAtomList(alMeshi, tor.atom4)); 
				int torsionNum = tor.getTorsionResNum();
				double weight = torRestrainSpringConst * (10.0*Math.PI/180.0) / pp[torsionNum][11];
				outString += (restrainTorString + ind1 + " " + ind2 + " " + ind3 + " " + ind4 + " " + (weight*Math.PI/180.0) + " " + 
						(pp[torsionNum][7]*180.0/Math.PI - 1e-4) + " " + (pp[torsionNum][7]*180.0/Math.PI + 1e-4) + " " + "\n");
			}
		}
		return outString;
	}
	
	public void addStringToKeyFile(String initialKeyFile, String outputKeyFile, String addString) {
		String[] lines = File2StringArray.f2a(initialKeyFile);
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outputKeyFile)));	
			for (int c=0 ; c<lines.length ; c++) {
				pw.println(lines[c]);
			}
			pw.println();
//			System.out.println(addString);
			pw.println(addString);
			pw.close();
		}
		catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * This method will create a new XYZ file where connectivity in chain breaks is (justifily) removed.
	 */
	public void fixBreaksInXYZ(String newXYZfile) {
		String[] lines = File2StringArray.f2a(tinkerXYZ);
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(newXYZfile)));	
			pw.println(lines[0]);
			for (int c=1 ; c<lines.length ; c++) {
				StringTokenizer st = new StringTokenizer(lines[c]);
				int tinkerInd = Integer.parseInt(st.nextToken());
				pw.print(tinkerInd + " ");
				pw.print(st.nextToken() + " "); //The atom type 
				pw.print(st.nextToken() + " "); //x
				pw.print(st.nextToken() + " "); //y 
				pw.print(st.nextToken() + " "); //z
				pw.print(st.nextToken() + " "); //the opls type
				while (st.hasMoreTokens()) {
					int connectInd = Integer.parseInt(st.nextToken());
					if ((tinker2meshi(tinkerInd)==-1) || (tinker2meshi(connectInd)==-1)) {
						pw.print(connectInd + " ");
					}
					else if (alMeshi.atomAt(tinker2meshi(tinkerInd)).chain().equals(alMeshi.atomAt(tinker2meshi(connectInd)).chain()) &&
							(Math.abs(alMeshi.atomAt(tinker2meshi(tinkerInd)).residueNumber() - alMeshi.atomAt(tinker2meshi(connectInd)).residueNumber())<2)) {
						pw.print(connectInd + " ");						
					}
					else {
//						System.out.println("Remove connectivity:\n" + 
//								alMeshi.atomAt(tinker2meshi(tinkerInd)) + "\n" +
//								alMeshi.atomAt(tinker2meshi(connectInd)));
					}
				}
				pw.println();
			}
			pw.close();
		}
		catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
	

	public static void main(String[] args) {
		String commandString = args[0];
		String refName = args[1];
		String doName = args[2];
		double disCutOff = Double.parseDouble(args[3]);
		String keyFile = doName+".key";
		MeshiProgram.initRandom(0);
		CommandList commands = new CommandList(commandString);
		DunbrackLib lib = new DunbrackLib(commands,0.99,100);
		TINKERandMESHIcomparison tinkRef = new TINKERandMESHIcomparison(refName+".xyz", refName+".pdb");
		tinkRef.setInactiveFarFromInterface(disCutOff);
//		for (int c =1000 ; c<2000 ; c++)
//			System.out.println(c + " " + tinkRef.meshi2tinker(c) + " " + tinkRef.alMeshi.atomAt(c));
//		System.exit(0);
		System.out.println("Doing: " + doName);			
		TINKERandMESHIcomparison tink = new TINKERandMESHIcomparison(doName+".xyz", doName+".pdb");
		tink.fixBreaksInXYZ(doName+".tink.xyz");
		tink.setInactiveFarFromInterface(tinkRef);
		tink.addStringToKeyFile(keyFile, doName+".tink.key", tink.getInactiveString());
		tink.addStringToKeyFile(doName+".tink.key", doName+".tink.key", tink.getRestrainString());
		MESHIonTriC meshiOnTriC  = new MESHIonTriC(doName+".pdb", lib );	
		double[][] pp = meshiOnTriC.getNearestRotInfo();
		tink.addStringToKeyFile(doName+".tink.key", doName+".tink.key", tink.getTorsionRestrainString(pp));
	}
	
}
