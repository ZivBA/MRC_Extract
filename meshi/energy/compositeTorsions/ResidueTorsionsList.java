package meshi.energy.compositeTorsions;

import java.util.Iterator;

import meshi.geometry.DistanceMatrix;
import meshi.geometry.TorsionList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.util.MeshiList;
import meshi.util.Updateable;
import meshi.util.UpdateableException;
import meshi.util.filters.Filter;

/** An updateable list of ResidueTorsions.
 * 
 * @author El-ad David Amir
 *
 */
public class ResidueTorsionsList
	extends MeshiList
	implements Updateable {

	/* updates counter for compatability with Updateable */
	int numberOfUpdates;
	
	public ResidueTorsionsList() {
		super( new IsResidueTorsions() );
	}
	
	/** Creates ResidueTorsionsList from protein. */
	public ResidueTorsionsList( Protein protein, DistanceMatrix dm ) {
		this();
		
		/* create list of all torsions in protein */
		TorsionList torsionList =
			TorsionList.createQuickAndDirtyTorsionList( protein, dm );
		
		Residue res;
		Iterator resIter = protein.residueIterator();
		
		/* create all residue torsions and add them to list */
		while( (res = (Residue) resIter.next()) != null ) {
			if( res.dummy() )
				continue;
			
			ResidueTorsions resTorsions = new ResidueTorsions( res, torsionList );
			add( resTorsions );
		}
		
		tagPreProline();
		
	}
	
	public void update(int numberOfUpdates) throws UpdateableException {
		if (numberOfUpdates == this.numberOfUpdates+1) {
			for ( int i = 0; i < size(); i++ )
				((ResidueTorsions) elementAt(i)).update( numberOfUpdates );
			this.numberOfUpdates++;
		}
		else if ( numberOfUpdates != this.numberOfUpdates ) 
			throw new RuntimeException("incorrect update number" );
	}

	
	/** A method to filter residues that have {PHI,PSI} torsions. 
	 * This method removes residues at chain termini or Ca traces, but keeps 
	 * incomplete residues.*/
	public ResidueTorsionsList filterPhiPsiResidues() {
		return (ResidueTorsionsList) filter(new PhiPsiResidueFilter(), new ResidueTorsionsList());
	}

    public static class PhiPsiResidueFilter implements Filter,CompositeTorsionsDefinitions {
	public boolean accept(Object obj) {
	    ResidueTorsions tors = (ResidueTorsions) obj;
	    return ((tors.getTorsion(PHI)!=null) && 
	    		(tors.getTorsion(PSI)!=null));
	}
    }

    
	/** A method to filter for incomplete residues. Note that the resulting list
	 * CONTAINS all the INCOMPLETE residues of the protein. It is inverse of the
	 * filterCompleteResidues() method.
	 **/
	public ResidueTorsionsList filterIncompleteResidues() {
		return (ResidueTorsionsList) filter(new IncompleteResidueFilter(), new ResidueTorsionsList());
	}
    
    public static class IncompleteResidueFilter implements Filter {
    	public boolean accept(Object obj) {
    	    ResidueTorsions tors = (ResidueTorsions) obj;
    	    return !tors.hasAllTorsions();
    	}
        }
    
	
	/** A method to filter residues that have all torsions */
	public ResidueTorsionsList filterCompleteResidues() {
		return (ResidueTorsionsList) filter(new CompleteResidueFilter(), new ResidueTorsionsList());
	}

    public static class CompleteResidueFilter implements Filter {
	public boolean accept(Object obj) {
	    ResidueTorsions tors = (ResidueTorsions) obj;
	    return tors.hasAllTorsions();
	}
    }

	/** Straightforward filter class to identify ResidueTorsions class. */
	private static class IsResidueTorsions implements Filter {
		public boolean accept(Object obj) {
			return (obj instanceof ResidueTorsions); 
		}		
	}

	/**
	 *A method to tag as Pre-Proline all the residues that are such.
	 */
	private void tagPreProline() {
		for (int c=0 ; c<size() ; c++) {
			ResidueTorsions thisRT = (ResidueTorsions) elementAt(c);
			int resnum = thisRT.getResidueNumber();
			ResidueTorsions nextRT = findResidueInList(resnum+1);
			if ((nextRT!=null) && (nextRT.getResidueType()==12) && //It's a pre-Pro 
				(thisRT.getResidueType()!=12) && (thisRT.getResidueType()!=5)) //This is not gly or pro
				thisRT.setToPP();
		}
	} 
	
	public ResidueTorsions findResidueInList(int num) {
		for (int c=0 ; c<size() ; c++) 
			if (((ResidueTorsions) elementAt(c)).getResidueNumber() == num)
				return (ResidueTorsions) elementAt(c);
		return null;
	}	
}
