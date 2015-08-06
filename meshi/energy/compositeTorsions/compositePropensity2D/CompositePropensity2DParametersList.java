package meshi.energy.compositeTorsions.compositePropensity2D;

import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.energy.compositeTorsions.SplinedPolynomialsLoader;
import meshi.parameters.Residues;

public class CompositePropensity2DParametersList extends ParametersList
	implements Residues,CompositeTorsionsDefinitions {
	
    private SplinedPolynomialsLoader spl = null;

	public CompositePropensity2DParametersList( String splinedPolynomialsFileName ) {
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
		for( int aac=0; aac<=OMNI; aac++ )
			add( createParameters( aac ) );		
	}
	
	public Parameters createParameters( int residueType ) {
		return new CompositePropensity2DParameters( residueType, spl );
	}

	/** Returns parameters for ResidueTorsions object according to
	 * residue type. */
	public Parameters parameters(Object Obj) {
		ResidueTorsions residueTorsions = (ResidueTorsions) Obj;
		
		/* scan internal array got possible residue types */
		for( int i=0; i<size(); i++ ) {
			CompositePropensity2DParameters cpp = 
				(CompositePropensity2DParameters) elementAt(i);
			if( cpp.getResidueType() == residueTorsions.getResidueType() )
				return cpp;
		}
		
		return null;
	}
	
	public Parameters createParameters(String line) {
		throw new RuntimeException( "CompositePropensity2DEnergy term uses " +
			"SplinedPolynomialsLoader in order to create parameters list" );
	}

}
