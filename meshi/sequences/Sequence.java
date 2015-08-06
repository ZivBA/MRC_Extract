package meshi.sequences;
import java.util.Iterator;

/**
 * Sequence is a single row alignment.
 **/
public class Sequence extends SequenceAlignment{
    public final SequenceCharFilter charFilter;
    public static final char UNKNOWN = '_';
    
        public Sequence(SequenceCharFilter charFilter) {
	super();
	this.charFilter = charFilter;
    }

    public Sequence(String comment, SequenceCharFilter charFilter) {
	this(charFilter);
	comments.add(comment);
    }

    public Sequence(String sequence, String comment, SequenceCharFilter charFilter) {
	this(charFilter);
	Character weirdChar = weirdChar(sequence);
	if (weirdChar == null) { // that is all characters in the sequence are either valid residue 
		                 // codes, gap or wildcard
		comments.add(comment);
		int number = 1;
		for (int i = 0; i < sequence.length(); i++) {
		    SequenceAlignmentCell cell;
	    	    char chr = sequence.charAt(i);
	    	    if (chr == SequenceAlignmentCell.GAP_CHAR) 
		 	cell = new SequenceAlignmentCell();
	    	    else {
		       cell = new SequenceAlignmentCell(chr,number);
		       number++;
	            }
	            SequenceAlignmentColumn column = new SequenceAlignmentColumn(cell);
	            add(column);
	        }
	}
    }
    
    
    public boolean add(Object obj) {
	SequenceAlignmentColumn column = (SequenceAlignmentColumn) obj;
	char c = column.getChar(0);
	Character cc = new Character(c);
	if (!charFilter.accept(cc))
	    throw new RuntimeException("Cannot add "+column.getChar(0)+" to "+getClass());
	return super.add(column);
    }

    public Character weirdChar(String sequence) {
	for (int i = 0; i < sequence.length(); i++) {
	    char c = sequence.charAt(i);
	    if ((c != SequenceAlignmentCell.GAP_CHAR) &
	        (c != SequenceAlignmentCell.WILDCARD_CHAR)) {
		Character cc = new Character(c);
		if (!charFilter.accept(cc)) return cc;
	    }
	}
	return null;
    }
    
    public String comment() { 
	if (comments.size() > 0) return (String) comments.elementAt(0);
	throw new RuntimeException("No comment to this Sequence object - "+
				   "probably an indication of a problem.");
    }


    public Sequence(Sequence source, SequenceCharFilter charFilter) {
	this(source.comment(),charFilter);
	add(source);
    } 	
    
    public int startsIn() {
	    int i;
	    boolean found =false;
	    for (i = 0;i <size() & (! found); i++)
		if (!cell(i).gap())found = true;
	    return (cell(i).number-i);
    }

    public String toString() {
	if (comments.size() == 0) {
	    if (size() == 0)
		throw new RuntimeException("Weird empty sequence without comment");
	    throw new RuntimeException("Weird sequence without comment");
	}
	int commentEnd = comment().indexOf("starts in ");
	if  (commentEnd == -1) commentEnd = comment().length();
	String out = "> "+comment().substring(0,commentEnd)+"; starts in "+startsIn();
	for (int i = 0; i < size(); i++) {
	    if(i%50 == 0) out += "\n";
	    out += charAt(i);
	}
	out+="*"+"\n";
	return out;
    }

    
    public SequenceAlignmentCell cell(int index) {
	if (index > size() -1) 
	    throw new RuntimeException ("No index "+index+
					" in Sequence instance of length "+size());
	AlignmentColumn column = columnAt(index);
	return (SequenceAlignmentCell) column.cell(0);
    }

    public char  getChar(int index) {
	return cell(index).getChar();
	
    }
    public String  getCharAsString(int index) {
	return cell(index).getCharAsString();
	
    }
    
    public char charAt(int index) {
	return cell(index).getChar();
    }

    public Sequence renumber(Sequence reference) {
	SequenceAlignment alignment = SequenceAlignment.identityAlignment(reference, this);
	SequenceAlignmentColumn firstColumn = (SequenceAlignmentColumn) alignment.elementAt(0);
	if (reference.size() == 0) return null;
	    if (! firstColumn.isExactMachWithGaps()) return renumber(reference.tail());
	if (! alignment.isExactMachWithGaps()) 
	    throw new RuntimeException("Cannot renumber \n"+this+"\n"+
				       "with\n"+
				       alignment);
	return renumber(alignment);
    }

    public Sequence tail() {
	Sequence out = new Sequence(comment(),charFilter);
	Iterator columns = iterator();
	if (! columns.hasNext()) return out;
	columns.next();
	while (columns.hasNext())
	    out.add(columns.next());
	return out;
    }
	
    public Sequence renumber(SequenceAlignment alignment) {
	AlignmentCell  pivot = null;
	int shift = 999999;
	for (int i = 0; (i < size()) & (pivot == null); i++) {
	    SequenceAlignmentCell  cell = cell(i);
	    for (Iterator columns = alignment.iterator(); (columns.hasNext()) & (pivot == null);) {
		SequenceAlignmentColumn  column = (SequenceAlignmentColumn) columns.next();
		if ((column.cell(1) == cell) & 
		    (!cell.gap())) {
		    pivot = column.cell(0);
		    shift = pivot.number - cell.number;
		}
	    }
	}
	if (pivot == null) 
	    throw new RuntimeException("Cannot renumber \n"+this+"\n"+
				       "with\n"+
				       alignment);

	Sequence out = new Sequence(comment(),charFilter);
	for (int i = 0; (i < size()) ; i++) {
	    SequenceAlignmentCell  cell = cell(i);
	    if (cell.gap()) out.add(cell);
	    else out.add(new SequenceAlignmentCell(cell.getChar(),cell.number+shift));
	}
	return out;
    }
    
    public boolean add(SequenceAlignmentCell cell) {
	return add(new SequenceAlignmentColumn(cell));
    }

}
