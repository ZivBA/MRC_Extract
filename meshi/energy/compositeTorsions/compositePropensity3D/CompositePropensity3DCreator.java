package meshi.energy.compositeTorsions.compositePropensity3D;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.compositeTorsions.ResidueTorsionsList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.parameters.MeshiPotential;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.KeyWords;

/** The propensity 3D energy calculates the proensity of each amino acid to be in a specific position 
 *  at the (phi,psi,chi1) torsion space. Note that if the residue being evaluated is GLY or ALA, or if
 *  its structure does not have chi1 defined (for example incomplete side chain), than its 2D propensity 
 *  is calculated instead. 
 */
public class CompositePropensity3DCreator
	extends EnergyCreator
	implements KeyWords, MeshiPotential, Residues {

	public CompositePropensity3DCreator(double weight) {
	super(weight);
    }
    
    public CompositePropensity3DCreator() {
		super( 1.0 );
	}
	
	public AbstractEnergy createEnergyTerm(Protein protein,
			DistanceMatrix distanceMatrix, CommandList commands) {
		/* retrieve parameters */
		if (parametersList == null) {
		    // Appanding the path to the filename list of the parameters.
	    	String[] strlist = new String[COMPOSITE_PROPENSITY_3D_PARAMETERS.length];
	    	String pathname = parametersDirectory(commands).concat("/");
	    	for (int cc=0 ; cc<COMPOSITE_PROPENSITY_3D_PARAMETERS.length ; cc++)
	        	strlist[cc] = pathname.concat(COMPOSITE_PROPENSITY_3D_PARAMETERS[cc]);
	        parametersList = new CompositePropensity3DParametersList(strlist);
	    }
		
		/* create residue torsions list for protein */
		ResidueTorsionsList rtl = (new ResidueTorsionsList(protein, distanceMatrix)).filterPhiPsiResidues();

		/* return energy */
		return new CompositePropensity3DEnergy(rtl, distanceMatrix,
		 (CompositePropensity3DParametersList) parametersList, weight(), "prop3D" );
	}	
}
