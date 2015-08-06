package meshi.energy.compositeTorsions.ramachandranAndChi1;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.compositeTorsions.ResidueTorsionsList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.parameters.MeshiPotential;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/** The Phi,Psi,Chi1 energy calculates the probability of each amino acid to be in a specific position 
 *  at the (phi,psi,chi1) torsion space. Note that if the residue being evaluated is GLY or ALA, or if
 *  its structure does not have chi1 defined (for example incomplete side chain), than the 2D probability 
 *  is calculated instead. 
 */
public class RamachandranAndChi1Creator
	extends EnergyCreator
	implements KeyWords, MeshiPotential, Residues {

	public RamachandranAndChi1Creator(double weight) {
	super(weight);
    }
    
    public RamachandranAndChi1Creator() {
		super( 1.0 );
	}
	
	public AbstractEnergy createEnergyTerm(Protein protein,
			DistanceMatrix distanceMatrix, CommandList commands) {
		/* retrieve parameters */
		if (parametersList == null) {
			String rsplFileName = 
				parametersDirectory(commands)+"/"+COMPOSITE_TORSIONS_PARAMETERS;
			parametersList = new RamachandranAndChi1ParametersList(rsplFileName);
		}
		
		/* create residue torsions list for protein */
		ResidueTorsionsList rtl = (new ResidueTorsionsList(protein, distanceMatrix)).filterPhiPsiResidues();

		/* return energy */
		return new RamachandranAndChi1Energy(rtl, distanceMatrix,
		 (RamachandranAndChi1ParametersList) parametersList, weight(), "Ramach3D" );
	}	
}
