package meshi.applications.TriC.yeast;

import meshi.applications.TriC.ReadTriCofOrganism;
import meshi.applications.TriC.TricAlignment;
import meshi.applications.TriC.TricYeastAlignment;
import meshi.applications.TriC.VectorOfOrganisms;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;

public class FindContacts  extends MeshiProgram implements Residues,AtomTypes {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		init(args);
		TricYeastAlignment alignmentYeast = new TricYeastAlignment();
		TricAlignment alignmentBov = new TricAlignment();
		VectorOfOrganisms vec = new VectorOfOrganisms();
		vec.add(new ReadTriCofOrganism());
		vec.add(new ReadTriCofOrganism("G:\\Home_Lehavim_28_11_2011\\TRiC\\Organisms\\human.txt"));
		vec.add(new ReadTriCofOrganism("G:\\Home_Lehavim_28_11_2011\\TRiC\\Organisms\\zebrafish.txt"));
		vec.add(new ReadTriCofOrganism("G:\\Home_Lehavim_28_11_2011\\TRiC\\Organisms\\ciona.txt"));
		vec.add(new ReadTriCofOrganism("G:\\Home_Lehavim_28_11_2011\\TRiC\\Organisms\\drosoph.txt"));
		vec.add(new ReadTriCofOrganism("G:\\Home_Lehavim_28_11_2011\\TRiC\\Organisms\\Celegance.txt"));
		vec.add(new ReadTriCofOrganism("G:\\Home_Lehavim_28_11_2011\\TRiC\\Organisms\\arabidopsis.txt"));
		vec.add(new ReadTriCofOrganism("G:\\Home_Lehavim_28_11_2011\\TRiC\\Organisms\\neospora.txt"));
		vec.add(new ReadTriCofOrganism("G:\\Home_Lehavim_28_11_2011\\TRiC\\Organisms\\candida.txt"));
		vec.add(new ReadTriCofOrganism("G:\\Home_Lehavim_28_11_2011\\TRiC\\Organisms\\mold.txt"));
		vec.add(new ReadTriCofOrganism("G:\\Home_Lehavim_28_11_2011\\TRiC\\Organisms\\plasmodium.txt"));
		vec.add(new ReadTriCofOrganism("G:\\Home_Lehavim_28_11_2011\\TRiC\\Organisms\\paramecium.txt"));
		vec.add(new ReadTriCofOrganism("G:\\Home_Lehavim_28_11_2011\\TRiC\\Organisms\\yeast.txt"));
		double[][][] SBs = new double[9][1000][1000];
		double[][][] HBs = new double[9][1000][1000];
		double[][][] BBs = new double[9][1000][1000];
		double[][][] CCs = new double[9][1000][1000];
		String letter1 = "AGZQHEBDI";
		String letter2 = "GZQHEBDAJ";
		String[] files = {"PAIR_AG_corrected_refined.pdb", 
				"PAIR_GZ_corrected_refined.pdb", 
				"PAIR_ZQ_corrected_refined.pdb", 
				"PAIR_QH_corrected_refined.pdb", 
				"PAIR_HE_corrected_refined.pdb", 
				"PAIR_EB_corrected_refined.pdb", 
				"PAIR_BD_corrected_refined.pdb", 
				"PAIR_DA_corrected_refined.pdb", 
				"1Q3R_refined.pdb"};
		for (int ind = 0; ind<9 ; ind++) {
			System.out.println("Reading: " + files[ind]);
			findContacts(files[ind], SBs[ind], HBs[ind], BBs[ind], CCs[ind]);
		}
		for (int resB=16 ; resB<1000 ; resB++) {
			boolean printed = false;
				for (int resA=16 ; resA<1000 ; resA++) {
				boolean printedUnits = false;
				for (int ind = 0; ind<9 ; ind++) {
					int resAunit;
					int resBunit;
					if (ind!=8) {
						resAunit = alignmentYeast.getNewResNum('K', resA, letter1.charAt(ind));
						resBunit = alignmentYeast.getNewResNum('K', resB, letter2.charAt(ind));
					}
					else {
						resAunit = resA;
						resBunit = resB;						
					}
					if ((resAunit!=-1) & (resBunit!=-1)) { // Both valid units
						String typeInter = "";
						if (SBs[ind][resAunit][resBunit]<999) {
							typeInter = "SB";
							printed=true;
							printedUnits=true;
						} else if (HBs[ind][resAunit][resBunit]<999) {
							printed=true;
							typeInter = "HH";
							printedUnits=true;
						} else if (CCs[ind][resAunit][resBunit]<999) {
							printed=true;
							typeInter = "CC";
							printedUnits=true;
						} else if (BBs[ind][resAunit][resBunit]<999) {
							printed=true;
							typeInter = "bb";				
							printedUnits=true;
						}  
						if (typeInter.length() > 0) {
							if (ind != 8) {
								System.out.println(resA + "  " + resAunit + "  " +  letter1.charAt(ind) + "  " + 
										"(" + alignmentYeast.getAAinSeq('K', resA) + ")" + 
										vec.getProfile(letter1.charAt(ind), alignmentBov.getNewResNum('K', resA, 'A')-1) + " " +
										" -> " +
										typeInter + 
										" <- " + 
										"(" + alignmentYeast.getAAinSeq('K', resB) + ")" +
										vec.getProfile(letter2.charAt(ind), alignmentBov.getNewResNum('K', resB, 'A')-1) + "  " +
										letter2.charAt(ind) + "  " + resBunit + "  " + resB);
							}
							else {
								System.out.println(resA + "  " + resAunit + "  " +  letter1.charAt(ind) + "  " + 
										"(" + alignmentYeast.getAAinSeq('K', resA) + ")" + 
										"=============" + " " +
										" -> " +
										typeInter + 
										" <- " + 
										"(" + alignmentYeast.getAAinSeq('K', resB) + ")" +
										"=============" + "  " +
										letter2.charAt(ind) + "  " + resBunit + "  " + resB);
							}
						}
					} // Both valid units
				} // Subunit type
				if (printedUnits) {
					System.out.println();
				}
			}
			if (printed) {
				System.out.println("------------------------------------------------------------");
			}
		}
		

//	char unitA = 'D';
//	char unitB = 'A';
//	
//		for (int c=0 ; c<1000 ; c++) {
//			for (int d=0 ; d<1000 ; d++) {
//				if (SBs[7][c][d]<999) {
//					System.out.println("SB: " + unitA + " " + c + " " + alignmentYeast.getNewResNum(unitA , c, 'K') + " " +
//							unitB + " " + d + " " + alignmentYeast.getNewResNum(unitB, d, 'K') + " " + SBs[7][c][d]);
//				}
//				else if (HBs[7][c][d]<999) {
//					System.out.println("HB: " + unitA + " " + c + " " + alignmentYeast.getNewResNum(unitA, c, 'K') + " " +
//							unitB + " " + d + " " + alignmentYeast.getNewResNum(unitB, d, 'K') + " " + HBs[7][c][d]);
//				}
//				else if (CCs[7][c][d]<999) {
//					System.out.println("CC: " + unitA + " " + c + " " + alignmentYeast.getNewResNum(unitA, c, 'K') + " " +
//							unitB + " " + d + " " + alignmentYeast.getNewResNum(unitB, d, 'K') + " " + CCs[7][c][d]);
//				}
//				else if (BBs[7][c][d]<999) {
//					System.out.println("BB: " + unitA + " " + c + " " + alignmentYeast.getNewResNum(unitA, c, 'K') + " " +
//							unitB + " " + d + " " + alignmentYeast.getNewResNum(unitB, d, 'K') + " " + BBs[7][c][d]);
//				}
//			}
//		}	
	}
	
	
	public static void findContacts(String pdbComplex, 
		double[][] SBs, double[][] HBs, double[][] BBs, double[][] CCs) {
		double SB_HB_th = 3.1;
		double C_C_th = 4.0;
		for (int c=0 ; c<1000 ; c++) {
			for (int d=0 ; d<1000 ; d++) {
				BBs[c][d] = SBs[c][d] = HBs[c][d] = CCs[c][d] = 1000.0;
			}
		}
		AtomList list = new AtomList(pdbComplex);
		for (int a1=0 ; a1<list.size() ; a1++) {
			for (int a2=a1+1 ; a2<list.size() ; a2++) {
				Atom atom1 = list.atomAt(a1);
				Atom atom2 = list.atomAt(a2);
				if (atom1.chain().charAt(0)!=atom2.chain().charAt(0)) {
					if ((!atom1.isBackbone | (atom1.residueName().equals("ALA") & atom1.name().equals("CB")) |
							(atom1.residueName().equals("GLY") & atom1.name().equals("CA"))) & 
						(!atom2.isBackbone | (atom2.residueName().equals("ALA") & atom2.name().equals("CB")) |
								(atom2.residueName().equals("GLY") & atom2.name().equals("CA")))) {
						// Check SBs and HBs
						if (((atom1.isNitrogen & atom2.isOxygen) |
								(atom2.isNitrogen & atom1.isOxygen)) |
							(atom1.isOxygen & atom2.isOxygen & 
							(atom1.residueName().equals("THR") |
							 atom1.residueName().equals("SER") |
							 atom1.residueName().equals("TYR") |
							 atom2.residueName().equals("THR") |
							 atom2.residueName().equals("SER") |
							 atom2.residueName().equals("TYR")))){
							// Check SB
							if (((atom1.residueName().equals("ARG") |
							atom1.residueName().equals("LYS") |
							atom1.residueName().equals("HIS")) & 
							(atom2.residueName().equals("ASP") |
							atom2.residueName().equals("GLU"))) |
							((atom2.residueName().equals("ARG") |
							atom2.residueName().equals("LYS") |
							atom2.residueName().equals("HIS")) & 
							(atom1.residueName().equals("ASP") |
							atom1.residueName().equals("GLU")))) {
								double dis = atom1.distanceFrom(atom2);
								if (dis < SB_HB_th) {
									if (SBs[atom1.residueNumber()][atom2.residueNumber()] > dis) {
										SBs[atom1.residueNumber()][atom2.residueNumber()] = dis;
									}
								}	
							}
							else { // Check HB
								double dis = atom1.distanceFrom(atom2);
								if (dis < SB_HB_th) {
									if (HBs[atom1.residueNumber()][atom2.residueNumber()] > dis) {
										HBs[atom1.residueNumber()][atom2.residueNumber()] = dis;
									}
								}	
							}
						}
						// Check CC
						if ((atom1.isCarbon | atom1.isSulfur) & (atom2.isCarbon | atom2.isSulfur)) {
							double dis = atom1.distanceFrom(atom2);
							if (dis < C_C_th) {
								if (CCs[atom1.residueNumber()][atom2.residueNumber()] > dis) {
									CCs[atom1.residueNumber()][atom2.residueNumber()] = dis;
								}
							}								
						}
					}
					//Check BB
					if ((atom1.isBackbone & !atom2.isBackbone) |
						(!atom1.isBackbone & atom2.isBackbone)) {
						if ((atom1.isNitrogen & atom2.isOxygen) |
							(atom2.isNitrogen & atom1.isOxygen)) {
							double dis = atom1.distanceFrom(atom2);
							if (dis < SB_HB_th) {
								if (BBs[atom1.residueNumber()][atom2.residueNumber()] > dis) {
									BBs[atom1.residueNumber()][atom2.residueNumber()] = dis;
								}
							}	
						}
					}
				}
			}
		}
	}
	
	
	protected static void init(String[] args) {
		int zvl = ALA; // force the reading of "meshi.parameters.Residues"
		zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"
		initRandom(333);
	}	
	

}
