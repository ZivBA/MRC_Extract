package meshi.energy.compositeTorsions;

import java.util.Iterator;

import meshi.geometry.Torsion;
import meshi.geometry.TorsionList;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Residue;
import meshi.parameters.Residues;
import meshi.util.Updateable;
import meshi.util.UpdateableException;

/** Allows easy access to all of a given residue's torsions.
 * 
 * @author El-ad David Amir
 *
 */
public class ResidueTorsions
	implements Updateable, CompositeTorsionsDefinitions, Residues {

	private Residue residue;
	private boolean isPP = false;
	
	/* the torsion angles array. notice that the array indices correspond
	 * to the torsion angle types, therefore, torsionAngles[0] is always
	 * null.
	 */
	private Torsion[] torsionAngles;
	
	/* updates counter for compatability with Updateable */
	private int numberOfUpdates = 0;
	
	/** stores all of the residue's torsion angles.
	 * 
	 * @param residue for which torsion angles are stored
	 * @param torsionList from where the torsion angles are taken
	 */
	public ResidueTorsions( Residue residue, TorsionList torsionList ) {
		if( residue.dummy() )
			throw new RuntimeException( "unable to process dummy residue" );
		
		this.residue = residue;
		
		/* create the torsion angles array */
		torsionAngles = new Torsion[TOTAL_TORSION_ANGLES];
		
		/* scan torsions list for this residue's torsions */
		Torsion tor;
		Iterator torIter = torsionList.iterator();
		while( (tor = (Torsion) torIter.next()) != null )
			/* check whether torsion belongs to this residue */
			if( tor.getTorsionResNum() == residue.number )
				/* check whether it's an "interesting" torsion angle */
				if( interestingTorsion( tor ) )
					/* and update array */
					torsionAngles[ torsionNum( tor ) ] = tor;
	}
	
	/** returns torsion of type torsionType. */
	public Torsion getTorsion( int torsionType ) {
		return torsionAngles[ torsionType ];
	}
	
	/** returns residue type */
	public int getResidueType() {
		return residue.type;
	}
	
	/** returns residue's number */
	public int getResidueNumber() {
		return residue.number;
	}
	
	/** returns residue */
	public Residue getResidue() {
		return residue;
	}
	
	public void setToPP() {isPP=true;}

	public boolean isPP() {return isPP;}
	
	/** verifies all the torsions used by this residue are available.
	 * This includes phi,psi torsion angles.
	 */
	public boolean hasAllTorsions() {
		/* all residues should have phi/psi torsion angles */
		if( torsionAngles[PHI] == null || torsionAngles[PSI] == null )
			return false;
		
		/* verify residues with at least one possible sidechain torsion */
		if( NUM_SIDECHAIN_TORSIONS[getResidueType()] >= 1 )
			if( torsionAngles[CHI_1] == null )
				return false;
		
		/* verify residues with at least two possible sidechain torsions */
		if( NUM_SIDECHAIN_TORSIONS[getResidueType()] >= 2 )
			if( torsionAngles[CHI_2] == null )
				return false;
				
		/* verify residues with at least three possible sidechain torsions */
		if( NUM_SIDECHAIN_TORSIONS[getResidueType()] >= 3 )
			if( torsionAngles[CHI_3] == null )
				return false;

		/* verify residues with at least four possible sidechain torsions */
		if( NUM_SIDECHAIN_TORSIONS[getResidueType()] >= 4 )
			if( torsionAngles[CHI_4] == null )
				return false;

		return true;
	}
	
	/** Applies force to torsion angles. Force is applied to each
	 * atom of the torsion.
	 */
	public void applyForce( int torsionType, double force ) {
		Torsion tor = getTorsion( torsionType );
		
		if( !tor.atom1.frozen() ) {
		    tor.atom1.addToFx(force * tor.dTorsionDx1());
		    tor.atom1.addToFy(force * tor.dTorsionDy1());
		    tor.atom1.addToFz(force * tor.dTorsionDz1());
		}
		
		if( !tor.atom2.frozen() ) {
		    tor.atom2.addToFx(force * tor.dTorsionDx2());
		    tor.atom2.addToFy(force * tor.dTorsionDy2());
		    tor.atom2.addToFz(force * tor.dTorsionDz2());
		}
		
		if( !tor.atom3.frozen() ) {
		    tor.atom3.addToFx(force * tor.dTorsionDx3());
		    tor.atom3.addToFy(force * tor.dTorsionDy3());
		    tor.atom3.addToFz(force * tor.dTorsionDz3());
		}
		
		if( !tor.atom4.frozen() ) {
		    tor.atom4.addToFx(force * tor.dTorsionDx4());
		    tor.atom4.addToFy(force * tor.dTorsionDy4());
		    tor.atom4.addToFz(force * tor.dTorsionDz4());
		}
	}
	
	/** returns a list of all atoms in this residue torsions */
	public AtomList getAtoms(int[] whichTorsions) {
		AtomList atoms = new AtomList();
		
		for( Torsion tor : torsionAngles )
			if( tor != null )  
				for (int findTor : whichTorsions)
					if (torsionNum(tor) == findTor) {
						atoms.add( tor.atom1 );
						atoms.add( tor.atom2 );
						atoms.add( tor.atom3 );
						atoms.add( tor.atom4 );
					}
		return atoms;
	}
	
	public void update(int numberOfUpdates) throws UpdateableException {
		if( numberOfUpdates == this.numberOfUpdates+1 ) {
			/* advance one update for all torsions */
			for( Torsion tor : torsionAngles )
				if( tor != null ) {
					tor.update( numberOfUpdates );
				}
			this.numberOfUpdates++;
		}
		else if( numberOfUpdates != this.numberOfUpdates )
				throw new RuntimeException( "incorrect update number" );
	}

	/** Converts attributes of ResidueTorsions to string. */
	public String toString() {
		/* convert torsions to string */
		String torsionAnglesStr = "";
		for( Torsion tor : torsionAngles )
			if( tor != null )
				torsionAnglesStr += " " + tor.torsion();
		
		return "residue: " + residue.name + "(" + residue.number +
				") torsionAngles: " + torsionAnglesStr;		
	}
	
	/** Partial conversion to string. Convert phi, psi and chi_1
	 * torsion angles. */
	public String toStringPhiPsiChi1() {
		String outStr = "";
		
		outStr += residue.name + "(" + residue.number + ") ";
		if( torsionAngles[PHI] == null || torsionAngles[PSI] == null )
			return outStr;
		
		outStr += torsionAngles[PHI].torsion() + " " + torsionAngles[PSI].torsion();
		
		if( torsionAngles[CHI_1] == null )
			return outStr;
		
		outStr += " " + torsionAngles[CHI_1].torsion();

		return outStr;
	}

	/** Partial conversion to string. Convert phi, psi
	 * torsion angles. */
	public String toStringPhiPsi() {
		String outStr = "";
		
		outStr += residue.number + " " + residue.type + " ";
		
		outStr += 180.0/Math.PI*torsionAngles[PHI].torsion() + " " + 180.0/Math.PI*torsionAngles[PSI].torsion();
		
		return outStr;
	}
	
	/** checks whether a given torsion angle is in [PHI,PSI,CHI_1-CHI_4]. */
	private static boolean interestingTorsion( Torsion tor ) {
		return ( torsionNum( tor ) != UNIDENTIFIED_TORSION_TYPE );
	}
	
	/** returns type of given torsion */
	private static int torsionNum( Torsion tor ) {
		if( tor == null )
			return UNIDENTIFIED_TORSION_TYPE;
		if( tor.name().equals( "PHI" ) )
			return PHI;
		else if( tor.name().equals( "PSI" ) )
			return PSI;
		else if( tor.name().equals( "CHI1" ) )
			return CHI_1;
		else if( tor.name().equals( "CHI2" ) )
			return CHI_2;
		else if( tor.name().equals( "CHI3" ) )
			return CHI_3;
		else if( tor.name().equals( "CHI4" ) )
			return CHI_4;
		else if( tor.name().equals( "OMG" ) )
			return OMG;
		else
			return UNIDENTIFIED_TORSION_TYPE;
	}
}
