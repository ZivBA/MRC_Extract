package programs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.StringTokenizer;

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


public class ExtractEVparameters extends MeshiProgram implements Residues, AtomTypes{ 
 
    public static void main(String[] args){
    	
    	// ----------------------------------------
    	// These are really for the user to define:
    	double cutoff = 4.5;
    	int howmany = 10000000;
    	// ----------------------------------------
    	
    	init(args); 
    	    	
    	int totalCounter = 0;
    	Protein model = null;
    	DistanceMatrix dm = null;
		int[] atomicTypeConverter = null;
		int maxAtomType=-1;
		int tmp;
		BufferedReader br;
		BufferedWriter bw;
		StringTokenizer stok;
		String line = "";
		int tsaiAtom1 = -1 ,tsaiAtom2 = -1;
	    DecimalFormat fmt = new DecimalFormat("0.##");    
		

    	// Converting the 190 MESHI atom types to the 14 mentioned in Tsai 99'
    	// -------------------------------------------------------------------
    	String Meshi2TsaiFileName = "C:/Users/Nir/Nir_Programs_eclipse/MESHI_E/meshicurrent/parameters/meshiPotential/SolvateMESHI2Tsai.dat";
   		System.out.println("Reading parameter file: " + Meshi2TsaiFileName);
    	try{
    		// first pass on the file - to find the maximal atom type
    		br = new BufferedReader(new FileReader(Meshi2TsaiFileName));
    	    line = br.readLine();
    	    while (line != null) {
    	    	stok = new StringTokenizer(line);
    	    	tmp = Atom.type(stok.nextToken().trim());
    	    	if (tmp>maxAtomType)
    	    	   maxAtomType = tmp;
    	    	line = br.readLine();
    	    }
    	    br.close();
    	    atomicTypeConverter = new int[maxAtomType+1];
    	    for(int c=0 ; c<atomicTypeConverter.length; c++)
    	    	atomicTypeConverter[c] = -1; 
    		// second pass on the file - reading the new types
    		br = new BufferedReader(new FileReader(Meshi2TsaiFileName));
    	    line = br.readLine();
    	    while (line != null) {
    	    	stok = new StringTokenizer(line);
    	    	tmp = Atom.type(stok.nextToken().trim());
    	    	atomicTypeConverter[tmp] = Integer.valueOf(stok.nextToken().trim()).intValue()-1;
    	    	line = br.readLine();
    	    }
    	    br.close();
    	}
    	catch(Exception e) {
    	    throw new RuntimeException(e.getMessage());
    	}
    	
    	// Going over the models
    	String[] models = File2StringArray.f2a("C:/Users/Nir/Loop_Building_Project/listPDBs.txt");
    	try{
    		bw = new BufferedWriter(new FileWriter("distances.txt"));
  
    	for (int i=0 ; (i<models.length) && (totalCounter<howmany) ; i++) { 		
			System.out.println("ReadingTrying to minimize: " + models[i]);
			model = new Protein(new AtomList(models[i]), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
			dm = new DistanceMatrix(model.atoms(), 5.5, 2.0, 4);
			DistanceList dl = dm.nonBondedList();
			for (int c=0 ; c<dl.size(); c++) {
				if (dl.distanceAt(c).distance()<cutoff) {
					tsaiAtom1 = atomicTypeConverter[dl.distanceAt(c).atom1().type];
					tsaiAtom2 = atomicTypeConverter[dl.distanceAt(c).atom2().type];
					if ((tsaiAtom1!=14)&&(tsaiAtom2!=14))  {  // These are NOT hydrogens, which have a Tsai ID of 14
						totalCounter++;
						bw.write(99 + " " + tsaiAtom1 + " " + tsaiAtom2 + " " + fmt.format(dl.distanceAt(c).distance()) + "\n");
						System.out.println(99 + " " + tsaiAtom1 + " " + tsaiAtom2 + " " + fmt.format(dl.distanceAt(c).distance()));
					}
				}	
			}
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
