package programs;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.DecimalFormat;

import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;


public class ExtractEVallTypes extends MeshiProgram implements Residues, AtomTypes{ 
 
    public static void main(String[] args){
    	
    	// ----------------------------------------
    	// The bins:
    	double[] bins = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5};
    	// ----------------------------------------
    	double[][][] data;
    	
    	init(args); 
    	
    	data = new double[Atom.numberOfTypes()][Atom.numberOfTypes()][bins.length];
    	double cutoff = bins[bins.length-1]+0.5*(bins[bins.length-1]-bins[bins.length-2]);
    	int ty1,ty2,ttemp;
    	int ind;
    	double dis;

    	
    	// Going over the models
    	String[] models = File2StringArray.f2a("C:/Users/Nir/Loop_Building_Project/listPISCES.txt");
    	for (int i=0 ; i<models.length ; i++) { 		
			System.out.println("Reading: " + models[i]);
			Protein model = new Protein(new AtomList("PISCES/"+models[i]), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
			DistanceMatrix dm = new DistanceMatrix(model.atoms(), bins[bins.length-1]+1.0 , 2.0, 4);
			DistanceList dl = dm.nonBondedList();
			for (int c=0 ; c<dl.size(); c++) {
				dis = dl.distanceAt(c).distance();
				if ((dis<cutoff) & (dl.distanceAt(c).atom1().residueNumber() != dl.distanceAt(c).atom2().residueNumber())) {
					ty1 = dl.distanceAt(c).atom1().type;
					ty2 = dl.distanceAt(c).atom2().type;
/*					if ((((ty1==PCD)  && dl.distanceAt(c).atom2().isBackbone && dl.distanceAt(c).atom2().isOxygen) || 
							((ty2==PCD)  && dl.distanceAt(c).atom1().isBackbone && dl.distanceAt(c).atom1().isOxygen)) &&
							(Math.abs(dl.distanceAt(c).atom1().residueNumber()-dl.distanceAt(c).atom2().residueNumber()) < 2)) {
						System.out.println("-----------------------------------------------------------------------------------------");
						System.out.println(dl.distanceAt(c).atom1() + "\n" + dl.distanceAt(c).atom2());
					}
*/					
					if (ty2>ty1) {
						ttemp = ty1;
						ty1 = ty2;
						ty2 = ttemp;
					}
					if (dis>bins[bins.length-1])
						data[ty1][ty2][bins.length-1]++;
					else {
						ind = 0;
						do {
							ind++;
						} while (dis>bins[ind]);
						data[ty1][ty2][ind] += (dis - bins[ind-1])/(bins[ind] - bins[ind-1]);
						data[ty1][ty2][ind-1] += (bins[ind] - dis)/(bins[ind] - bins[ind-1]);
					}
				}
			}
    	}
    	
    	// Outputting
    	try{
    	    DecimalFormat fmt = new DecimalFormat("0.##");    
    		BufferedWriter bw = new BufferedWriter(new FileWriter("distances.txt"));
    		for (int c1=0; c1<Atom.numberOfTypes() ; c1++)
        		for (int c2=0; c2<Atom.numberOfTypes() ; c2++) {
        			bw.write(c1 + " " + c2 + "   ");
        			System.out.print(c1 + " " + c2 + "   ");
        			for (int d=0 ; d<bins.length ; d++) {
        				bw.write(fmt.format(data[c1][c2][d]) + " ");
        				System.out.print(fmt.format(data[c1][c2][d]) + " ");
        			}
        			bw.write("\n");
        			System.out.println();
        		}
    	    bw.close();
    	}
    	catch(Exception e) {
    	    throw new RuntimeException(e.getMessage());
    	}
    }

    
    
     
    protected static void init(String[] args) {
 
	/**** NOTE *** the next two lines. Because of a BUG in the Java VM, the 
	 * interfaces "Residues" and "AtomTypes" are not loaded automatically when MinimizeProtein initialize. 
	 * For this purpose these two lines are crucial wherever these two interfaces are implemented. The user might 
	 * rightfully feel that these two lines are "black magic" programming, but happily to our knowledge this is 
	 * the only bizarre phenomenon we are aware of in meshi.
	 **/
	int zvl = ALA; // force the reading of "meshi.parameters.Residues"
	zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"

	initRandom(333);
    }
}
