package meshi.applications.Valpuesta;

import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;

public class AddAllAtomsToUNK extends MeshiProgram implements Residues,AtomTypes  {

	public static void main(String[] args) {
		init(args);

		AtomList all = new AtomList("2XSM_allALA.pdb");
		String outfile = "2XSM_allALA_full";
		String chains="ABCDEFGHIJKLMNOP";
		for (int let=0 ; let<chains.length() ; let++) {
			try {
				BufferedWriter bw = new BufferedWriter(new FileWriter(outfile+"_"+chains.charAt(let)+".pdb"));
				Protein full = new Protein(all.chainFilter(""+chains.charAt(let)), new ResidueExtendedAtoms(ADD_ATOMS));
				full.atoms().setChain(""+chains.charAt(let));
				for (int c=0 ; c<full.atoms().size() ; c++) {
					bw.write(full.atoms().atomAt(c).toString() + "\n");
				}
				bw.write("TER\n");			
				bw.write("END\n");			
				bw.close();
			}
			catch(Exception e) {
				throw new RuntimeException(e.getMessage());
			}
		}
			
}
	
		
	protected static void init(String[] args) {
		int zvl = ALA; // force the reading of "meshi.parameters.Residues"
		zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"
		initRandom(333);
	}	

	
	
}
