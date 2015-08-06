package meshi.applications.prediction;

import meshi.sequences.FastaList;
import meshi.sequences.ResidueSequence;
import meshi.sequences.Sequence;
import meshi.util.CommandList;
import meshi.util.Key;

public class PredictionResidueSequence extends ResidueSequence {
    private static Sequence temp = null;
    public PredictionResidueSequence(CommandList commands, Key key) {
	super();
	String fileName=commands.firstWord(key).secondWord();
	FastaList fastaList = new FastaList(fileName);
	Sequence temp = new ResidueSequence(fastaList);
	comments.add(temp.comment());
	add(temp);
    }
}
