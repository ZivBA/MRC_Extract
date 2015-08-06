package meshi.energy.distanceConstrains;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.DummyResidue;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.filters.Filter;

public class TemplateDistanceConstrainsCreator extends EnergyCreator implements KeyWords {
    private double intraSegmentFactor;
    private double intraSegmentTolerance;
    private double interSegmentFactor;
    private double interSegmentTolerance;
    private double saturation;
    private double unsatisfiedCutoff;
    private double upToCutoff;
    private int[] fragments;
    private String distanceConstrainsMask = null;
    private boolean debug = true;
    private static Protein protein;
    private Filter filter;

    public TemplateDistanceConstrainsCreator(Filter filter)  {
	super(TEMPLATE_DISTANCE_CONSTRAINS);
	this.filter = filter;
	if (filter == null) throw new RuntimeException("filter is null");
    }

    public  AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					   CommandList commands) {
	return new DistanceConstrainsEnergy(protein, (DistanceConstrainParametersList) parametersList, weight());
    }
	
	
	/**
	 *This method is used to filterout unsatisfied constrains. 
	 *Its output is a new DistanceConstrainParametersList that keeps only the satisfied 
	 *constrains. Also the inner parameterList of the Creator is updated.
	 **/
    public DistanceConstrainParametersList filterParametersList(Protein protein) {
	parametersList = new DistanceConstrainParametersList((DistanceConstrainParametersList) parametersList, 
							     protein, unsatisfiedCutoff);
    return (DistanceConstrainParametersList) parametersList;
    }

    public void setParametersList(Protein protein, CommandList commands) {
	this.protein = protein;

	CommandList constrainCommands = commands.firstWordFilter(key);

	intraSegmentFactor = constrainCommands.secondWord(INTRA_SEGMENT_FACTOR).thirdWordDouble();
	if (debug) System.out.println("intraSegmentFactor = "+intraSegmentFactor); 

	intraSegmentTolerance = constrainCommands.secondWord(INTRA_SEGMENT_TOLERANCE).thirdWordDouble();
	if (debug) System.out.println("intraSegmentTolerance = "+intraSegmentTolerance); 

	interSegmentFactor = constrainCommands.secondWord(INTER_SEGMENT_FACTOR).thirdWordDouble();
	if (debug) System.out.println("interSegmentFactor = "+interSegmentFactor); 

	interSegmentTolerance = constrainCommands.secondWord(INTER_SEGMENT_TOLERANCE).thirdWordDouble();
	if (debug) System.out.println("interSegmentTolerance = "+interSegmentTolerance); 

	saturation = constrainCommands.secondWord(SATURATION).thirdWordDouble();
	if (debug) System.out.println("saturation = "+saturation); 

	unsatisfiedCutoff = constrainCommands.secondWord(UNSATISFIED_CUTTOF).thirdWordDouble();
	if (debug) System.out.println("unsatisfiedCutoff = "+unsatisfiedCutoff); 

	upToCutoff = constrainCommands.secondWord(UP_TO_CUTOFF).thirdWordDouble();
	if (debug) System.out.println("upToCutoff = "+upToCutoff); 

	distanceConstrainsMask = getDistanceConstrainsMask(commands);
	System.out.println("# No constraints are taken from residues marked with X: \n"+ 
			   distanceConstrainsMask.substring(1)); //To "hide" the first dummy residue (number 0)

	//weightWasSet = true;

	System.out.println(getFragments(protein));
	parametersList = new DistanceConstrainParametersList(protein, fragments, 
							     intraSegmentFactor,
							     intraSegmentTolerance,
							     interSegmentFactor,
							     interSegmentTolerance,
							     saturation,
							     upToCutoff,
							     distanceConstrainsMask,
							     filter);
    }

    // -------------- get fragments ---------------
    private String getFragments(Protein protein) {
	System.out.println(" getting fragments of "+protein);
	fragments = new int[protein.residues().size()];
	int fndr = protein.residues().firstNonDummyResidueNumber();
	System.out.println("First Non-dummy residue "+
			   " "+protein.residues().elementAt(fndr));
	int fragment = 0;
	for (int ires = 0; ires < fndr; ires++)
	    fragments[ires] = -1;
	fragments[fndr] = 0;
	String out = "\nFragments:\n"+protein.residues().elementAt(fndr)+"\t"+fragments[fndr]+"\n";
	for (int ires = fndr+1; ires < protein.residues().size(); ires++) {
	    if (protein.residues().elementAt(ires) instanceof DummyResidue) 
		fragments[ires] = -1;
	    else {
		if ((((Residue)protein.residues().elementAt(ires-1)).ca() == null) |
		    (((Residue)protein.residues().elementAt(ires)).ca() == null)) {
		    fragment++;
		    fragments[ires] = fragment;
		}
		else {
		    Residue residueIm1 = (Residue) protein.residues().elementAt(ires-1);
		    Residue residueI = (Residue) protein.residues().elementAt(ires);
		    if ((residueIm1.ca() != null) & (residueI.ca() != null)) {
			double dis =  residueIm1.ca().distance(residueI.ca()).distance();
			if ((dis < 2.8) | (dis > 3.865)) { // the lower bound corresponds to cis peptide
			    fragment++;
			    fragments[ires] = fragment;
			}
			else fragments[ires] = fragment;
		    }
		}
	    }
	    out += protein.residues().elementAt(ires)+"\t"+fragments[ires]+"\n";
	}
	return out;
    }



    public static String getDistanceConstrainsMask(CommandList commands) {
	String dub = "";
	String fragment = "";
	String line;
	int cc;
	int firstNonDummyResidueIndex = protein.residues().firstNonDummyResidueIndex();
	for (cc=0 ; cc<firstNonDummyResidueIndex ; cc++) {
	    dub += "X";
	}
	if (commands.keyExists(DISTANCE_CONSTRAINS_MASK))
	    dub += commands.getSequence(DISTANCE_CONSTRAINS_MASK);
	else {
	    for (cc=firstNonDummyResidueIndex ; cc<protein.residues().size() ; cc++) {
		dub += "-";
	    }
	}
	    
	return dub;
    }
}
