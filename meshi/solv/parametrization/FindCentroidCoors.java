package meshi.solv.parametrization;

import meshi.geometry.Angle;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.Torsion;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;
import meshi.util.overlap.Overlap;

/**
 * 
 * The centroids are calculated, and the output is the code that goes into the 'buildCentroid'
 * method in 'ResidueBuilder'.
 **/

public class FindCentroidCoors extends MeshiProgram  implements Residues{ 

	// We get these user defined parameters from the command line
	// ----------------------------------------------------------
	// The structure database:
	static String listOfStructures = "listPDBs.txt";
//	static String listOfStructures = "C:/Users/Nir/Loop_Building_Project/listPISCES.txt";

	public static void main(String[] args){
		
		initRandom(333);


		int[] nHeavyAtom = {
				1, // ALA
				2, // CYS
				4, // ASP
				5, // GLU
				7, // PHE
				0, // GLY
				6, // HIS
				4, // ILE
				5, // LYS
				4, // LEU
				4, // MET
				4, // ASN
				3, // PRO
				5, // GLN
				7, // ARG
				2, // SER
				3, // THR
				3, // VAL
				10,// TRP
				8  // TYR
				};
		
		int[] constantIndices = {0,1,2};

		double[][][] shiftingCoors = new double[20][3][];
		for (int c=0 ; c<nHeavyAtom.length ; c++) {
			shiftingCoors[c][0] = new double[nHeavyAtom[c]+3];
			shiftingCoors[c][1] = new double[nHeavyAtom[c]+3];
			shiftingCoors[c][2] = new double[nHeavyAtom[c]+3];
		}
		
		double[][] NCACcoors = {{-1.459, 0, 0.546}, {0, 0, 1.424} , {0, 0, 0}};
		
		double[][] averageCoors = new double[20][3];
		double[] atomsNumber = new double[20];
		
						
		for (int doType = 0; doType<20; doType++) 
		if (doType!=5) {
			
		for (int c=0 ; c<20 ; c++)
			averageCoors[c][0] = averageCoors[c][1] = averageCoors[c][2] = atomsNumber[c] = 0.0; 
		
		Atom atomN,atomCA,atomC;
		Residue residue;
		String[] models = File2StringArray.f2a(listOfStructures);
		for (int i=0 ; i<models.length ; i++) { 		
			System.out.println("Reading: " + models[i]);
			Protein model = new Protein(new AtomList("proteins/" +models[i]), new ResidueExtendedAtoms(ADD_ATOMS));
			for (int res=0 ; res<model.residues().size(); res++) {
				if ((model.residues().residueAt(res).type == doType) &&
					(model.residues().residueAt(res).atoms().filter(new AtomList.NonHydrogen()).size() == (nHeavyAtom[doType]+4))) {
					residue = model.residues().residueAt(res);
					atomN = residue.atoms().findAtomInList("N", residue.number); 
					atomCA = residue.atoms().findAtomInList("CA", residue.number); 
					atomC = residue.atoms().findAtomInList("C", residue.number); 
					shiftingCoors[doType][0][0] = atomN.x();
					shiftingCoors[doType][1][0] = atomN.y();
					shiftingCoors[doType][2][0] = atomN.z();
					shiftingCoors[doType][0][1] = atomCA.x();
					shiftingCoors[doType][1][1] = atomCA.y();
					shiftingCoors[doType][2][1] = atomCA.z();
					shiftingCoors[doType][0][2] = atomC.x();
					shiftingCoors[doType][1][2] = atomC.y();
					shiftingCoors[doType][2][2] = atomC.z();
					int counter = 3;
					for (int c=0 ; c<residue.atoms().size(); c++) {
						if (!residue.atoms().atomAt(c).isHydrogen && (!residue.atoms().atomAt(c).isBackbone || residue.atoms().atomAt(c).name().equals("CB"))) {
							shiftingCoors[doType][0][counter] = residue.atoms().atomAt(c).x();
							shiftingCoors[doType][1][counter] = residue.atoms().atomAt(c).y();
							shiftingCoors[doType][2][counter] = residue.atoms().atomAt(c).z();
							counter++;
						}
					}
					Overlap.rmsPartialAltRMS(NCACcoors, shiftingCoors[doType], constantIndices);
					for (int c=3 ; c<shiftingCoors[doType][0].length; c++) {
						averageCoors[doType][0] += shiftingCoors[doType][0][c];
						averageCoors[doType][1] += shiftingCoors[doType][1][c];
						averageCoors[doType][2] += shiftingCoors[doType][2][c];
						atomsNumber[doType]++;
					}
				}
			}
		}
		
		
		averageCoors[doType][0] = averageCoors[doType][0]/atomsNumber[doType];
		averageCoors[doType][1] = averageCoors[doType][1]/atomsNumber[doType];
		averageCoors[doType][2] = averageCoors[doType][2]/atomsNumber[doType];
		Protein model = new Protein(new AtomList("proteins/" +models[1]), new ResidueExtendedAtoms(ADD_ATOMS));
		for (int res=0 ; res<model.residues().size(); res++) {
			if (model.residues().residueAt(res).type == doType) {
				residue = model.residues().residueAt(res);
				residue.atoms().findAtomInList("N", residue.number).setXYZ(-1.459, 0, 0);
				residue.atoms().findAtomInList("CA", residue.number).setXYZ(0, 0, 0);
				residue.atoms().findAtomInList("C", residue.number).setXYZ(0.546, 1.424, 0);
				residue.atoms().findAtomInList("CB", residue.number).setXYZ(averageCoors[doType][0],
						averageCoors[doType][1],
						averageCoors[doType][2]);
				DistanceMatrix dm = new DistanceMatrix(model.atoms(), 5.5 , 1.0, 4);
				Angle angle2 = new Angle(residue.atoms().findAtomInList("C", residue.number),
						residue.atoms().findAtomInList("CA", residue.number),
						residue.atoms().findAtomInList("CB", residue.number),
						dm,
						false);
				Angle angle1 = new Angle(residue.atoms().findAtomInList("N", residue.number),
						residue.atoms().findAtomInList("C", residue.number),
						residue.atoms().findAtomInList("CA", residue.number),
						dm,
						false);
				double dis = residue.atoms().findAtomInList("CA", residue.number).distanceFrom(residue.atoms().findAtomInList("CB", residue.number));
				Torsion tor = new Torsion(angle1,angle2,dm);
				System.out.println("   case "+doType+":  bond = "+dis+"; //"+Residue.nameOneLetter(doType)+"\n"+
						"            angle = "+angle2.angle()+";\n"+
						"            tor = " + tor.torsion()+";\n"+
				"            break;");
				break;
			}				
		}		
		
	} // Do type loop

// This will show where the centroids are on a protein:
// ----------------------------------------------------		
//		String[] models = File2StringArray.f2a(listOfStructures);
//		Protein model = new Protein(new AtomList(models[0]), new ResidueExtendedAtoms(ADD_ATOMS));
//		for (int res=0 ; res<model.residues().size(); res++) {
//			if (model.residues().residueAt(res).type<20)
//			ResidueBuilder.buildCentroid(model.residues().residueAt(res));
//		}
//		try {
//			model.atoms().print(new MeshiWriter(models[0]+".centroid.pdb"));
//		}
//		catch (Exception e) {
//			System.out.print("\nThere was a problem writing:\n" + e + "\n\nContinuing...\n\n");
//		}

	} // Of main
		
}




	
