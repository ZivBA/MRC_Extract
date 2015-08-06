package meshi.energy.compositeTorsions.smoothRotamerLibrary;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.compositeTorsions.ResidueTorsionsList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.parameters.MeshiPotential;
import meshi.util.CommandList;
import meshi.util.KeyWords;

public class SmoothRotamerLibraryEnergyCreator
	extends EnergyCreator
	implements KeyWords, MeshiPotential {

	public SmoothRotamerLibraryEnergyCreator(double weight) {
	super(weight);
    }
    
    public SmoothRotamerLibraryEnergyCreator() {
		super( 1.0 );
	}
	
	public AbstractEnergy createEnergyTerm(Protein protein,
			DistanceMatrix distanceMatrix, CommandList commands) {
		/* retrieve parameters */
		String srlplFileName = 
			parametersDirectory(commands)+"/"+COMPOSITE_TORSIONS_PARAMETERS;
		SmoothRotamerLibraryParametersList srlpl = 
			new SmoothRotamerLibraryParametersList(srlplFileName);
		
		/* create residue torsions list for protein */
		ResidueTorsionsList rtl = (new ResidueTorsionsList(protein, distanceMatrix )).filterCompleteResidues(); 		
		
		/* return energy */
		return new SmoothRotamerLibraryEnergy(
				rtl, distanceMatrix, srlpl, weight(), "smooth" );
	}

}
