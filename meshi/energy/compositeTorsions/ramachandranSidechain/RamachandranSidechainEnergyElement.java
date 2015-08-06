package meshi.energy.compositeTorsions.ramachandranSidechain;

import meshi.energy.EnergyElement;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.parameters.Residues;

/** Encapsulation of RamachandranSidechain energy value for a single residue.
 * Much like the parameters for this energy function, each residue type
 * has its own class of energy element.
 * @author El-ad David Amir
 *
 */
public abstract class RamachandranSidechainEnergyElement
	extends EnergyElement
	implements Residues, CompositeTorsionsDefinitions {

	protected ResidueTorsions residueTorsions;
	protected RamachandranSidechainParameters rsp;
	protected double weight;

	public RamachandranSidechainEnergyElement(
			ResidueTorsions residueTorsions,
			RamachandranSidechainParameters rsp,
			double weight ) {
		this.residueTorsions = residueTorsions;
		this.rsp = rsp;
		this.weight = weight;
		
		if( !legalResidueType() )
			throw new RuntimeException( "energy element residue type mismatch" );
		
		setAtoms();
		updateFrozen();
	}

	protected void setAtoms() {
		int[] interestingTorsions = {PHI,PSI,CHI_1,CHI_2,CHI_3,CHI_4};
		atoms = residueTorsions.getAtoms(interestingTorsions);		
	}

	/** Reports energy values. Currently switched off. */
	protected void monitor( double energy, double ... derivs ) {
		if( false ) {
			System.out.println( "[[BEGIN]] RamachandranSidechainEnergy monitor information" );
			System.out.println( residueTorsions );
			System.out.println( "total energy = " + energy );
			for( int i=0; i<derivs.length; i++ )
				System.out.println( "deriv #" + (i+1) + " = " + derivs[i] );
			System.out.println( "[[END]] RamachandranSidechainEnergy monitor information" );
		}
	}
	
	/** verifies residue type is a legal residue types for class. */
	protected abstract boolean legalResidueType();
	
}
