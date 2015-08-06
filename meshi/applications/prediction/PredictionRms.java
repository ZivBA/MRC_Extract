package meshi.applications.prediction;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.sequences.AtomAlignment;
import meshi.sequences.ResidueAlignment;
import meshi.util.CommandList;
import meshi.util.Rms;
import meshi.util.filters.Filter;
import meshi.util.filters.KolDichfin;

public class PredictionRms extends Rms {
    private static boolean debug = true;
    public PredictionRms(ResidueAlignment residueAlignmet) {
	super(residueAlignmet);
    }
	
    public PredictionRms(AtomAlignment atomAlignment) {
	super(atomAlignment);
    }

    public static double rms(Protein protein0, CommandList commands) {
	return rms(protein0, commands, new KolDichfin()); // KolDichfin - a filter that accept all
    }

    public static double rms(Protein protein0, CommandList commands, Filter filter) {
	AtomAlignment atomAlignment = getAtomAlignment(protein0, commands, filter);
	if (atomAlignment == null) return -1;
	Rms rms = new Rms(atomAlignment);
	return rms.getRms();
    }

    public static double rms(Protein protein0, Protein protein1) {
    	ResidueAlignment residueAlignment = new ResidueAlignment(protein0, protein1);
	AtomAlignment atomAlignment = new AtomAlignment(residueAlignment, new KolDichfin());
	if (atomAlignment == null) return -1;
	Rms rms = new Rms(atomAlignment);
	return rms.getRms();
    }

    public static double gdt(Protein protein0, CommandList commands) {
	return gdt(protein0, commands, new KolDichfin()); // KolDichfin - a filter that accept all
    }

    public static double gdt(Protein model, CommandList commands, Filter filter) {
	try {
		CommandList alignmentCommands = commands.firstWordFilter(SUPERIMPOSE);
		String referenceFileName = alignmentCommands.secondWord(REFERENCE).thirdWord();
		if (referenceFileName.equals(NONE)) {
			 	System.out.println("No reference structure");
				return -1;
		}
                Protein reference = new Protein(referenceFileName);
		int refLength = reference.atoms().CAFilter().size();
		ResidueAlignment residueAlignment = new ResidueAlignment(reference,model);
	        AtomAlignment atomAlignment = new AtomAlignment(residueAlignment,filter);
		PredictionAtomList referenceAtoms = new PredictionAtomList(atomAlignment, 0);
		PredictionAtomList modelAtoms = new PredictionAtomList(atomAlignment, 1);

	        AtomList newReferenceAtoms = referenceAtoms.duplicate();
	        newReferenceAtoms.renumber();
	        AtomList newModelAtoms = modelAtoms.duplicate();
	        newModelAtoms.renumber();
	        if (true)
	        	throw new RuntimeException("\n\nNir disables this function 5.4.2006\n\n");
	        return GDTcalculator.gdt(newReferenceAtoms, newModelAtoms);//, refLength);
	}
	catch (Exception ex) {System.out.println("\n\n"+"Failed to calculate gdt\n"+ex);}
	return -1;
    }
	    


    private static AtomAlignment getAtomAlignment(Protein protein0, CommandList commands, Filter filter) {
	CommandList alignmentCommands = commands.firstWordFilter(SUPERIMPOSE);
	String referenceFileName = alignmentCommands.secondWord(REFERENCE).thirdWord();
	if (referenceFileName.equals(NONE)) {
	    System.out.println("No reference structure");
	    return null;
	}
	else {
	    Protein protein1 = new Protein(referenceFileName);
	    if (debug) System.out.println("protein1 = "+protein1);
	    String mode = alignmentCommands.secondWord(MODE).thirdWord();
	    if (!mode.equals(ALL_CA.key)) 
		throw new RuntimeException("Currently the only implemented mode of "+SUPERIMPOSE+"\n"+
					   "is "+ALL_CA);
	    System.out.println("Generating ResidueAlignment of \n"+protein1.sequence()+"\n"+protein0.sequence());
	    ResidueAlignment residueAlignment = new ResidueAlignment(protein1, protein0);
	    return  new AtomAlignment(residueAlignment, filter);
	}
    }
}
