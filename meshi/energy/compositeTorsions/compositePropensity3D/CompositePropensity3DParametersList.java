package meshi.energy.compositeTorsions.compositePropensity3D;

import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.energy.compositeTorsions.SplinedPolynomialsLoader;
import meshi.parameters.Residues;

public class CompositePropensity3DParametersList extends ParametersList
	implements Residues,CompositeTorsionsDefinitions {
	
    private SplinedPolynomialsLoader spl2D = null;
    private SplinedPolynomialsLoader spl3D = null;

	public CompositePropensity3DParametersList( String splinedPolynomialsFileName[]) {
		super();

		/* attempt to load polynomials from file */
		try {
		    System.out.println("Loading "+this+" parameters from "+splinedPolynomialsFileName[0]);
			spl2D = new SplinedPolynomialsLoader( splinedPolynomialsFileName[0] );
		    System.out.println("Loading "+this+" parameters from "+splinedPolynomialsFileName[1]);
			spl3D = new SplinedPolynomialsLoader( splinedPolynomialsFileName[1] );
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
		return new CompositePropensity2DParameters( residueType, spl2D );
	}
	public Parameters createParameters3D( int residueType ) {
		return new CompositePropensity3DParameters( residueType, spl3D );
	}

	/** Returns parameters for ResidueTorsions object according to
	 * residue type. */
	public Parameters parameters(Object Obj) {
		ResidueTorsions residueTorsions = (ResidueTorsions) Obj;
		
		if (residueTorsions.getTorsion( CHI_1 ) != null) {  // Giving a 3D parameters
			for( int i=0; i<size(); i++ ) {
				if ((elementAt(i) instanceof CompositePropensity3DParameters) && 
				( ((CompositePropensity3DParameters) elementAt(i)).getResidueType() ==
				 residueTorsions.getResidueType() ))
				 	return (Parameters) elementAt(i);
			}	
		}
		else { // Giving a 2D parameters
			for( int i=0; i<size(); i++ ) {
				if ((elementAt(i) instanceof CompositePropensity2DParameters) && 
				( ((CompositePropensity2DParameters) elementAt(i)).getResidueType() ==
				 residueTorsions.getResidueType() ))
				 	return (Parameters) elementAt(i);
			}				
		}
		return null;
	}
	
	public Parameters createParameters(String line) {
		throw new RuntimeException( "CompositePropensity3DParameterList term uses " +
			"SplinedPolynomialsLoader in order to create parameters list" );
	}

}
