package meshi.applications.corpus;

import meshi.energy.EnergyCreator;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.util.CommandList;

public class CorpusQuickNoSolvRot1 extends Corpus {

	public CorpusQuickNoSolvRot1(String PDBfile, CommandList commands,
			EnergyCreator[] energyCreators) {
		super(PDBfile, commands, energyCreators);
		// TODO Auto-generated constructor stub
	}

	public CorpusQuickNoSolvRot1(String exsitingCorpusFile) {
		super(exsitingCorpusFile);
		// TODO Auto-generated constructor stub
	}

	
	// Overriding this slow method so that all the SolvRot1 energies are 0.
	// ************************
	protected void calculateSolRot1(String PDBfile , CommandList commands) {
		// making the auxilary protein, and putting it in Rot1
		AtomList al1 = new AtomList(PDBfile);
		al1.moveCMtoOrigin();
		al1.renumber();
		Protein aux = new Protein(al1 , new ResidueExtendedAtoms(ADD_ATOMS));	
	    aux.defrost();
	    
		for (int c=1; c<(Nres-1) ; c++)
		if ((resNum[c-1]==resNum[c]-1) && (resNum[c+1]==resNum[c]+1)) { // This residue is not at the termini
			int num = resNum[c];
			System.out.println("AAAAAA: " + num);
			for (int mutateTo=0 ; mutateTo<20 ; mutateTo++)
				energies[c][1][mutateTo] = 1e-100;
		}
	}
	
	
}
