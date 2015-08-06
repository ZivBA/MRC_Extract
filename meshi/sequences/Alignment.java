package meshi.sequences;
import java.util.Iterator;

import meshi.util.MeshiList;
import meshi.util.filters.Filter;
import meshi.util.string.StringList;

/**
 * An alignment of proteins. Various sub-classes align different views of proteins (as sequences
 * lists of residues and lists of atoms.
 * Conceptually, an alignment is a matrix, where the rows are proteins and the columns align the 
 * corresponding elements of these proteins. The implementation however, is a bit different. 
 * It is implemented as a list of columns thus enjoying the rich functionality of MeshiList.
 **/
public class Alignment extends MeshiList{
    /**
     * A comment is associated with each protein.
     **/
    public final StringList comments;
    private Filter filter;
    /**
     * The filter determine the type of the alignment (sequence, residue or atom).
     **/
    public Alignment(Filter filter) {
	super(filter);
	this.filter = filter;
	comments = new StringList(); 
    }
    
    public AlignmentColumn columnAt(int index) {
	if (index > size()-1) throw new RuntimeException("No column number "+index+
						       " in Alignment instance of length "+
						       size());
	return (AlignmentColumn) elementAt(index);
    } 

    
    /**
     * Fetches a column with a specific cell number in a specific row. This column now is
     * a handle to the corresponding elements of the other proteins.
     **/
    public AlignmentColumn getColumn(int row, int number) { 
	for (Iterator columns = iterator(); columns.hasNext();) { 
	    AlignmentColumn column = (AlignmentColumn) columns.next(); 
	    if (column.cell(row).number == number) return column; 
	} 
	return null; 
    } 


    /**
     * Does the alignemnt include gaps? .
     **/
    public boolean hasGaps() { 
	for (Iterator iter = iterator(); iter.hasNext();) { 
	    AlignmentColumn column = (AlignmentColumn) iter.next(); 
	    if (column.hasGap()) return true; 
	} 
	return false; 
    } 
    
}
	
