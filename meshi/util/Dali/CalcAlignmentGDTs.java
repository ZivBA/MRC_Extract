package meshi.util.Dali;

import java.io.IOException;

import meshi.applications.prediction.GDTcalculator;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.MinimizerException;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;

public class CalcAlignmentGDTs implements Residues {

	public static void main(String[] args) throws MinimizerException, LineSearchException, IOException {
		String path = "/home/nirka/projects/directRef/"; // if not empty string must ends with '/'

		String protID = args[0].trim();
		MeshiProgram.initRandom(777);
		Protein nat = new Protein(path + "Natives/"+protID+".pdb" , new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		double gdt_temp = GDTcalculator.gdt(nat.atoms(), 
				new AtomList(path + "PreModels/"+protID+"_homo.pdb"));
		double gdt_full_nr = GDTcalculator.gdt(nat.atoms(), 
				new AtomList(path + "PreModels/"+protID+"_full.pdb"));
		double gdt_full = GDTcalculator.gdt(nat.atoms(), 
				new AtomList(path + "PreModels/"+protID+"_ref.pdb"));
		System.out.println("GDT TEMP: " + protID + "  " + gdt_temp);
		System.out.println("GDT FULL: " + protID + "  " + gdt_full_nr);
		System.out.println("GDT REFI: " + protID + "  " + gdt_full);

	}	
	
}
