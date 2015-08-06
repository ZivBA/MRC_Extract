package meshi.energy.compositeTorsions.ramachandranAndChi1;

import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.energy.compositeTorsions.SplinedPolynomialsLoader;
import meshi.parameters.Residues;

public class RamachandranAndChi1ParametersList extends ParametersList
	implements Residues,CompositeTorsionsDefinitions {
	
    private SplinedPolynomialsLoader spl = null;

	public RamachandranAndChi1ParametersList( String splinedPolynomialsFileName) {
		super();

		/* attempt to load polynomials from file */
		try {
		    System.out.println("Loading "+this+" parameters from "+splinedPolynomialsFileName);
			spl = new SplinedPolynomialsLoader( splinedPolynomialsFileName );
		}
		catch( Exception e ) {
			e.printStackTrace();
			throw new RuntimeException( "unable to read polynomials parameters file" );
		}
		
		/* create parameters for each amino acid type and OMNI */
		for( int aac=0; aac<=OMNI; aac++ ) {
			add( createParameters2D( aac ) );		
			add( createParameters3D( aac ) );
		}		
	}
	
	public Parameters createParameters2D( int residueType ) {
		return new RamachandranAndChi12DParameters( residueType, spl );
	}
	public Parameters createParameters3D( int residueType ) {
		return new RamachandranAndChi13DParameters( residueType, spl );
	}

	/** Returns parameters for ResidueTorsions object according to
	 * residue type. */
	public Parameters parameters(Object Obj) {
		ResidueTorsions residueTorsions = (ResidueTorsions) Obj;
		
		if (residueTorsions.getTorsion( CHI_1 ) != null) {  // Giving a 3D parameters
			for( int i=0; i<size(); i++ ) {
				if ((elementAt(i) instanceof RamachandranAndChi13DParameters) && 
				( ((RamachandranAndChi13DParameters) elementAt(i)).getResidueType() ==
				 residueTorsions.getResidueType() ))
				 	return (Parameters) elementAt(i);
			}	
		}
		else { // Giving a 2D parameters
			for( int i=0; i<size(); i++ ) {
				if ((elementAt(i) instanceof RamachandranAndChi12DParameters) && 
				( ((RamachandranAndChi12DParameters) elementAt(i)).getResidueType() ==
				 residueTorsions.getResidueType() ))
				 	return (Parameters) elementAt(i);
			}				
		}
		return null;
	}
	
	public Parameters createParameters(String line) {
		throw new RuntimeException( "RamachandranAndChi1ParameterList term uses " +
			"SplinedPolynomialsLoader in order to create parameters list" );
	}

}
