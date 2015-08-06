package programs;

import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.dssp.DSSP;
import meshi.util.file.File2StringArray;
import meshi.util.rotamericTools.RotamericTools;


public class ExtractPropensityData extends MeshiProgram implements Residues, AtomTypes{ 
 
    public static void main(String[] args){
    	  	
    	init(args); 
    	  	
    	// Going over the models
    	try{
    		BufferedWriter bw = new BufferedWriter(new FileWriter("phipsiData.txt"));
    		String[] prots = File2StringArray.f2a("C:/Users/Nir/Loop_Building_Project/extractingPropensity/pdbList.txt");
    		for (int i=0 ; i<prots.length ; i++) { 		
    			System.out.println("Reading: " + prots[i]);
    			Protein model = new Protein(new AtomList("../Pisces_REDUCE/"+prots[i]+".pdb"), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
    			DSSP dssp = new DSSP("Pisces_DSSP/"+prots[i]+".dssp");
    			double[][] pp = RotamericTools.phipsi(model, new DistanceMatrix(model.atoms(), 5.5, 2.0 , 4));
    			for (int c=0 ; c<model.residues().size() ; c++) {
    				int resNum = model.residueAt(c).number;
    				if (model.residueAt(c)!=null) {
    					if (model.residueAt(c).ca()!=null) {
    						if ((model.residue(resNum-1)!=null) && (model.residue(resNum+1)!=null) && (pp[resNum]!=null)) {
    							if (dssp.SSofRes(resNum,' ')=='H') { 
    								System.out.println(i + " " + resNum + " " + model.residueAt(c).type + " 2 " +
    										(pp[resNum][0]*180/Math.PI) + " " + (pp[resNum][1]*180/Math.PI) + " " + dssp.relACCofRes(resNum,' '));
    								bw.write(i + " " + resNum + " " + model.residueAt(c).type + " 2 " +
    										(pp[resNum][0]*180/Math.PI) + " " + (pp[resNum][1]*180/Math.PI) + " " + dssp.relACCofRes(resNum,' ') + "\n");
    							}
    							else if (dssp.SSofRes(resNum,' ')=='E') { 
    								System.out.println(i + " " + resNum + " " + model.residueAt(c).type + " 1 " +
    										(pp[resNum][0]*180/Math.PI) + " " + (pp[resNum][1]*180/Math.PI) + " " + dssp.relACCofRes(resNum,' '));
    								bw.write(i + " " + resNum + " " + model.residueAt(c).type + " 1 " +
    										(pp[resNum][0]*180/Math.PI) + " " + (pp[resNum][1]*180/Math.PI) + " " + dssp.relACCofRes(resNum,' ') + "\n");
    							}
    							else {
    								System.out.println(i + " " + resNum + " " + model.residueAt(c).type + " 0 " +
    										(pp[resNum][0]*180/Math.PI) + " " + (pp[resNum][1]*180/Math.PI) + " " + dssp.relACCofRes(resNum,' '));
    								bw.write(i + " " + resNum + " " + model.residueAt(c).type + " 0 " +
    										(pp[resNum][0]*180/Math.PI) + " " + (pp[resNum][1]*180/Math.PI) + " " + dssp.relACCofRes(resNum,' ') + "\n");
    							}								
    						}
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
