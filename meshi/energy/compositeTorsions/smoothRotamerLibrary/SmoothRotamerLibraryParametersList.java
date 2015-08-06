package meshi.energy.compositeTorsions.smoothRotamerLibrary;

import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.energy.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.compositeTorsions.ResidueTorsions;
import meshi.energy.compositeTorsions.SplinedPolynomialsLoader;
import meshi.parameters.Residues;

public class SmoothRotamerLibraryParametersList
	extends ParametersList
	implements CompositeTorsionsDefinitions, Residues {

	/** Loader contains all of the polynomials in the data file. */
	private static SplinedPolynomialsLoader spl = null;
	
	/** Create a new parameters file from the polynomials parameters data file. */
	public SmoothRotamerLibraryParametersList(
			String splinedPolynomialsFileName ) {
		super();

		/* do not load spl if it has already been loaded */
		if( spl == null )
			/* attempt to load polynomials from file */
			try {
				spl = new SplinedPolynomialsLoader( splinedPolynomialsFileName );
			}
			catch( Exception e ) {
				e.printStackTrace();
				throw new RuntimeException( "unable to read polynomials parameters file" );
			}
		
		/* create parameters for each amino acid type */
		for( int aac=0; aac<20; aac++ )
			add( createParameters( aac ) );
	}
	
	/** Returns parameters for ResidueTorsions object according to
	 * residue type. */
	public Parameters parameters(Object Obj) {
		ResidueTorsions residueTorsions = (ResidueTorsions) Obj;
		
		/* scan internal array got possible residue types */
		for( int i=0; i<size(); i++ ) {
			SmoothRotamerLibraryParameters srlp = 
				(SmoothRotamerLibraryParameters) elementAt(i);
			if( srlp.getResidueType() == residueTorsions.getResidueType() )
				return srlp;
		}
		
		return null;
	}

	/** Creates a SmoothRotamerLibraryParameters for residue type. */
	public Parameters createParameters( int residueType ) {
		/* classify residue according to type and call appropriate
		 * constructor class. types are grouped according to number
		 * of side chain torsion angles.
		 */
		if( NUM_SIDECHAIN_TORSIONS[residueType] == 0 )
			return new SmoothRotamerLibraryParametersChi0( residueType, spl );
		else if( NUM_SIDECHAIN_TORSIONS[residueType] == 1 )
			return new SmoothRotamerLibraryParametersChi1( residueType, spl );
		else if( NUM_SIDECHAIN_TORSIONS[residueType] == 2 )
			return new SmoothRotamerLibraryParametersChi2( residueType, spl );
		else if( NUM_SIDECHAIN_TORSIONS[residueType] == 3 )
			return new SmoothRotamerLibraryParametersChi3( residueType, spl );
		else if( NUM_SIDECHAIN_TORSIONS[residueType] == 4 )
			return new SmoothRotamerLibraryParametersChi4( residueType, spl );
		else
			throw new RuntimeException( "unrecognized amino acid type" );
	}

	/** Obsolete in SmoothRotamerLibraryParametersList. */
	public Parameters createParameters(String line) {
		throw new RuntimeException( "SmoothRotamerLibraryEnergy term uses " +
				"SplinedPolynomialsLoader in order to create parameters list" );
	}
}
