package meshi.applications.rotamerSearch;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

import programs.PutHydrogens;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.external.ComplexMESHIconversion;
import meshi.util.file.MeshiWriter;

public class SearchInterface extends MeshiProgram {
	
	public final double maxDisCBtoAtomInAnotherChain = 5.0; // All residues whose CB atoms are closer than this to another chain will be refined.
	public final int N_ITERATIONS = 6;
	public final boolean WITH_PERMUTE = true;
	
	private MultipleResidueSearch search;
	private Protein prot;

	public SearchInterface(CommandList commands , Protein prot , String writeInterfaceFile) {
		this.prot = prot;
		boolean[] interfaceResidues = interfacingResidues(prot , writeInterfaceFile);
		Vector<Residue> residues = new Vector<Residue>();
		for (int c=0 ; c<interfaceResidues.length ; c++) {
			if (interfaceResidues[c]) {
				residues.add(prot.residue(c));
			}
		}
		search = new MultipleResidueSearch(prot, commands, residues);
	}

	public void scan() {
		search.searchSequential(N_ITERATIONS, WITH_PERMUTE, randomNumberGenerator());
	}
	
	public void writeProt(String fileName, String chainA, String chainB) {
		Atom.resetNumberOfAtoms();
		AtomList list = ComplexMESHIconversion.MEHSI2complex(prot).duplicate();
		list.chainFilter("A").setChain("#");
		list.chainFilter("B").setChain("^");
		list.chainFilter("#").setChain(chainA);
		list.chainFilter("^").setChain(chainB);
		try {
			list.noOXTFilter().filter(new AtomList.NonHydrogen()).print(new MeshiWriter(fileName));
		} catch (IOException e) {
			throw new RuntimeException("Could not write file");
		}
	}
	
	
	/**
	 * The output is a boolean[2000] array (I assume that each subunit is not more 999 residues). True is interfacing.
	 */
	public boolean[] interfacingResidues(Protein workingProt, String writeInterfaceFile) {
		boolean[] out = new boolean[2000];
		for (int c=0 ; c<2000 ; c++) {
			out[c] = false;
		}
		// Looking from A to B
		for (int c=0 ; c<1000 ; c++) {
			Atom atom = workingProt.atoms().findAtomInList("CB", c);
			if ((atom!=null) && (!atom.residueName().equals("ALA"))
					&& (!atom.residueName().equals("PHE"))
					&& (!atom.residueName().equals("LEU"))
					&& (!atom.residueName().equals("ILE"))
					&& (!atom.residueName().equals("MET"))
					&& (!atom.residueName().equals("PRO"))
					&& (!atom.residueName().equals("VAL"))
					&& (!atom.residueName().equals("TRP"))
					&& (!atom.residueName().equals("TYR"))) {
				for (int d=0 ; d<workingProt.atoms().size() ; d++) {
					if (workingProt.atoms().atomAt(d).residueNumber()>1000) {
						if (atom.distanceFrom(workingProt.atoms().atomAt(d)) < maxDisCBtoAtomInAnotherChain) {
							out[c] = true;
						}
					}
				}
			}
		}
		// Looking from B to A
		for (int c=1000 ; c<2000 ; c++) {
			Atom atom = workingProt.atoms().findAtomInList("CB", c);
			if ((atom!=null) && (!atom.residueName().equals("ALA"))					
					&& (!atom.residueName().equals("PHE"))
					&& (!atom.residueName().equals("LEU"))
					&& (!atom.residueName().equals("ILE"))
					&& (!atom.residueName().equals("MET"))
					&& (!atom.residueName().equals("PRO"))
					&& (!atom.residueName().equals("VAL"))
					&& (!atom.residueName().equals("TRP"))
					&& (!atom.residueName().equals("TYR"))) {
				for (int d=0 ; d<workingProt.atoms().size() ; d++) {
					if (workingProt.atoms().atomAt(d).residueNumber()<1000) {
						if (atom.distanceFrom(workingProt.atoms().atomAt(d)) < maxDisCBtoAtomInAnotherChain) {
							out[c] = true;
						}
					}
				}
			}
		}
		int counter = 0;
		for (int c=0 ; c<2000 ; c++) {
			if (out[c])
				counter++;
		}
		System.out.println(counter + " residues in the interface.");
		try {
			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(writeInterfaceFile)));	
			for (int c=0 ; c<2000 ; c++) {
				if (out[c])
					pw.println(c);
			}
			pw.close();
		}
		catch (Exception e) {
			throw new RuntimeException(e);
		}		
		return out;
	}

	
	/**
	 * Will create an "A","B" complex of the desired pair.
	 */
	public static AtomList getPair(String fullComplexFileName , String chainA, String chainB) {
		AtomList fullList = new AtomList(fullComplexFileName);
		Atom.resetNumberOfAtoms();
		AtomList tmp1 = fullList.chainFilter(chainA).duplicate();
		tmp1.setChain("A");
		AtomList tmp2 = fullList.chainFilter(chainB).duplicate();
		tmp2.setChain("B");
		AtomList out = new AtomList();
		out.add(tmp1);
		out.add(tmp2);
		return out;
	}
	
	/**
	 * Will create a ready for work protein from an "A","B" atom list.
	 */
	public static Protein makeWorkableProtein(AtomList ABlist, CommandList commands) {
		Protein workingProt = ComplexMESHIconversion.complex2meshi(ABlist);
		PutHydrogens.adjustHydrogens(commands, workingProt);		
		return workingProt;
	}
	
	
	/**
	 * @param args
	 */		
	public static void main(String[] args) {
		initRandom(999);	
		CommandList commands = new CommandList(args[0]);
		Protein workProt = makeWorkableProtein(getPair(args[1], args[2], args[3]), commands);
		SearchInterface search = new SearchInterface(commands, workProt, args[4]);
		search.scan();
		search.writeProt(args[5], args[2], args[3]);
	}

	
}
