package meshi.applications.prediction;
import java.util.Iterator;

import meshi.molecularElements.Protein;
import meshi.sequences.Sequence;
import meshi.sequences.SequenceAlignment;
import meshi.sequences.SequenceAlignmentCell;
import meshi.sequences.SequenceAlignmentColumn;
import meshi.sequences.SequenceList;
import meshi.util.CommandList;
import meshi.util.Key;

public class PredictionSequenceAlignment extends SequenceAlignment {
    
    /**
     * A sequenceAlignment is extracted from a file whose name is referred to by the command line starting with key.
     * The first sequence of the alignment is assumed to be a subsequence of the reference sequence. The reference sequence is 
     * used to correct the start field of the aligned sequences. The reference itself is assummed to have correct start field 
     * (that is, it indicates the correct residue numbers). 
     **/
    private boolean debug = false;
    public PredictionSequenceAlignment(CommandList commands,Key key, Sequence reference) { 
	String fileName = commands.firstWord(key).secondWord();
	SequenceList sequenceList = new SequenceList(fileName);
	Sequence sequence0 = (Sequence) sequenceList.elementAt(0);
	Sequence sequence1 = (Sequence) sequenceList.elementAt(1);
	sequence0 = sequence0.renumber(reference);
	sequence1 = sequence1.renumber(new SequenceAlignment(sequence0,sequence1));
	SequenceAlignment temp = new SequenceAlignment(sequence0, sequence1);
	add(temp); /* This alignment is empty and it swallows all the elements (columns) of temp.*/ 
	comments.add(temp.comments); 
    } 
	
	
    public PredictionSequenceAlignment(Sequence targetSequence, Protein template, SequenceList sequenceList) { 
	super();
	Sequence templateSequence = template.sequence();
	Sequence targetFromAlignment = (Sequence) sequenceList.elementAt(0);
	Sequence templateFromAlignment = (Sequence) sequenceList.elementAt(1);
	SequenceAlignment templateAlignment = SequenceAlignment.identityAlignment(templateSequence, templateFromAlignment);
	if (! templateAlignment.isExactMachWithGaps()) throw new TemplateAlignmentException(template, sequenceList, templateAlignment);

	templateFromAlignment = templateFromAlignment.renumber(templateAlignment);

	SequenceAlignment tempAlignment = new SequenceAlignment(templateFromAlignment, targetFromAlignment);
	for (Iterator columns = tempAlignment.iterator(); columns.hasNext();) {
	    SequenceAlignmentColumn  column = (SequenceAlignmentColumn) columns.next();
	    SequenceAlignmentCell cell = (SequenceAlignmentCell) column.cell(1);
	    cell.addAttribute(column);
	}

	SequenceAlignment targetAlignment = SequenceAlignment.identityAlignment(targetSequence, targetFromAlignment);

	for (Iterator columns = targetAlignment.iterator(); columns.hasNext();) {
	    SequenceAlignmentColumn  column = (SequenceAlignmentColumn) columns.next();
	    SequenceAlignmentColumn  newColumn;
	    char gap = SequenceAlignmentCell.GAP_CHAR;
	    char cTarget = ((SequenceAlignmentCell)column.cell(0)).getChar();
	    char cSequence1 = ((SequenceAlignmentCell)column.cell(1)).getChar();
	    newColumn = new SequenceAlignmentColumn(2);
	    newColumn.add(0,column.cell(0));
	    if ((cSequence1 == gap) & (cTarget != gap)) 
		newColumn.add(1,new SequenceAlignmentCell()) ;
	    else {
		SequenceAlignmentCell tempCell = (SequenceAlignmentCell) column.cell(1);
		SequenceAlignmentColumn tempColumn = (SequenceAlignmentColumn) tempCell.getAttribute();
		newColumn.add(1,tempColumn.cell(0));
	    }
	    add(newColumn);
	}
	comments.add(templateSequence.comment());
	comments.add(targetSequence.comment());
    }

    private static class TemplateAlignmentException extends RuntimeException {
	public TemplateAlignmentException(Protein protein, SequenceList sequenceList, SequenceAlignment templateAlignment) {
	    super("\n"+"The template protein "+protein+" and the template sequence in the alignment file "+sequenceList.fileName()+"\n"+
		  "do not mach"+"\n"+templateAlignment);
	}
    }

}  
	/*
	--------------------------------------------------------------------------------------------------------------- 
	 This code was written originally by Chen, to allow for flexiblity in alignment formats. In the flexibility
	the template sequence could be either the first or the second, and the class would choose the correct one. This 
	feature was disabled by Nir on 1.5.2006 because it caused very strange bugs when the model was also the template.
	--------------------------------------------------------------------------------------------------------------- 
	if (templateAlignment0.score() > templateAlignment1.score()) {
	    templateAlignment = templateAlignment0;
	    targetFromAlignment = sequence1;
	    templateFromAlignment = sequence0;
	}
	else {
	    templateAlignment = templateAlignment1;
	    targetFromAlignment = sequence0;
	    templateFromAlignment = sequence1;
	}*/