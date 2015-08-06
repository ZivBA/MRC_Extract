package meshi.sequences;
import java.util.Iterator;
import java.util.LinkedList;

import meshi.util.CommandList;
import meshi.util.Key;
import meshi.util.KeyWords;
import meshi.util.MeshiList;
import meshi.util.filters.Filter;

public class SequenceList extends MeshiList implements KeyWords {
    private String fileName = null;
    protected static final SequenceCharFilter acceptAllChars = new AcceptAll(); 
    public SequenceList() {
	super(new IsSequenceFilter());
    }

    public SequenceList(String  fileName) {
	if (! FastaList.isFasta(fileName)) throw new RuntimeException("currently only Fasta files are used");
        SequenceList temp = new FastaList(fileName);	
	add(temp);
	this.fileName = fileName;
    }

    public SequenceList(SequenceAlignment alignment) {
	super(new IsSequenceFilter());
	SequenceAlignmentColumn firstColumn = (SequenceAlignmentColumn) alignment.elementAt(0); 
	int numberOfRows = firstColumn.size();
	for (int i = 0; i < numberOfRows; i++) {
	    Sequence sequence = new Sequence((String) alignment.comments.elementAt(i), 
					     acceptAllChars);
	    for (Iterator columns = alignment.iterator(); columns.hasNext();) {
		SequenceAlignmentColumn column = (SequenceAlignmentColumn) columns.next();
		sequence.add((SequenceAlignmentCell) column.cell(i));
	    }
	    add(sequence);
	}
    }

    //------------------------ extract sequences from the list ----------------------------
    private Sequence getSequence(String name) {
	    for (Iterator sequences = iterator(); sequences.hasNext();) {
		    Sequence sequence = (Sequence) sequences.next();
		    if (sequence.comment().indexOf(name) > -1) return sequence;
	    }
	    return null;
    }

    public ResidueSequence getResidueSequence(String name) {
	Sequence sequence = getSequence(name); 
	if (sequence == null) return null;
	return new ResidueSequence(sequence);
    }
    
    public ResidueSequence getResidueSequence() {
  	return getResidueSequence(AA_SEQUENCE.key);
    }
    
    public ResidueSequence getResidueSequence(CommandList commands, Key seqNameKey) {
	String name = commands.firstWord(seqNameKey).secondWord();
	return  getResidueSequence(name);
    }

    
    private static class IsSequenceFilter implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof Sequence);
	}
    }

    public void add(Sequence sequence) {
	    add(sequence, false); // do not flatten
    }

    public String toString() {
	    String out = "Sequence list \n";
	    for (Iterator sequences = iterator(); sequences.hasNext();) {
		    out += sequences.next().toString()+"\n";
	    }
	    return out;
    }

    public String fileName() { return fileName;}
    
    
    private static class AcceptAll extends SequenceCharFilter {
	public boolean accept(Object obj) { 
	    Character c = ((Character) obj).charValue();
	    return true;
	}
    }
    
    
}
	
