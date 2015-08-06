package meshi.energy.distanceConstrains;
import java.util.Iterator;

import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.molecularElements.Atom;
import meshi.molecularElements.DummyResidue;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.util.MeshiProgram;
import meshi.util.filters.Filter;

public class DistanceConstrainParametersList extends ParametersList{
    
    public DistanceConstrainParametersList(String parametersFileName) {
	super(parametersFileName,false);
    }

    public DistanceConstrainParametersList(DistanceConstrainParametersList original,
					    Protein protein,
					    double cutoff) {
	this(original, new SatisfiedConstrainsFilter(protein,cutoff));
    }

    public DistanceConstrainParametersList(DistanceConstrainParametersList original,
					    Filter filter) {
	super();
	System.out.println("Number of constrains before filtering: "+original.size());
	Iterator parametersIter = original.iterator();
	DistanceConstrainParameters parameters;
	while ((parameters = (DistanceConstrainParameters) parametersIter.next()) != null)
	    if (filter.accept(parameters)) add(parameters);
       System.out.println("Number of constrains after "+filter+": "+size());
    }

    public DistanceConstrainParametersList(Protein protein, int[] fragments, 
					   double intraSegmentFactor,
					   double intraSegmentTolerance,
					   double interSegmentFactor,
					   double interSegmentTolerance,
					   double distanceConsrainsSaturation,
					   double upToCutoff,
					   String distanceConstrainsMask,
					   Filter filter){
	super();
	double weight, tolerance;
	for (int ires = protein.residues().firstNonDummyResidueNumber(); ires < protein.residues().size(); ires++) {
	    Residue residueI = (Residue) protein.residues().elementAt(ires);
// 	    if (!(residueI instanceof DummyResidue)) System.out.println("xxxxx "+
// 									residueI.ca()+" "+filter.accept(residueI.ca()));
	    if ((!(residueI instanceof DummyResidue)) &&
		(residueI.ca() != null) &&
		filter.accept(residueI.ca()))
		for (int jres = ires+2; jres < protein.residues().size(); jres++) {
		    Residue residueJ = (Residue) protein.residues().elementAt(jres);
		    if ((!(residueJ instanceof DummyResidue))  &&
			(residueJ.ca() != null) &&
			filter.accept(residueJ.ca()))
			// The condition below is intended to reduce the number of constrains.
			// A more elegant way of doing this is desirable.
			if (((jres - ires < 10) & (fragments[ires] == fragments[jres])) |
			    (MeshiProgram.randomNumberGenerator().nextDouble() < distanceConsrainsSaturation)) {
			    if (fragments[ires] == fragments[jres]) {
				weight = intraSegmentFactor;
				tolerance = intraSegmentTolerance;
			    }
			    else { 
				weight = interSegmentFactor;
				tolerance = interSegmentTolerance;
			    }
			    Atom CaI = residueI.ca();
			    Atom CaJ = residueJ.ca();
			    if ((CaI != null) & (CaJ != null)) {
				double dis =  CaI.distanceFrom(CaJ);
				try {
				    if ((dis < upToCutoff) && 
					(distanceConstrainsMask.charAt(residueI.number) == '-') && 
					(distanceConstrainsMask.charAt(residueJ.number) == '-'))
					add(new DistanceConstrainParameters(residueI.name,
									    residueI.number,
									    "CA",
									    residueJ.name,
									    residueJ.number,
									    "CA",
									    dis, dis*tolerance, weight));
				}
				catch (StringIndexOutOfBoundsException e) { 
				    throw new RuntimeException("problem with template based "+
							       "DistanceConstrainParamerterList\n"+
							       "Apaprently distanceConstrainsMask is too short "+
							       "("+distanceConstrainsMask.length()+")\n"+
							       "while trying to add parameters for "+residueI+
							       " and "+residueJ+"\n"+e);
				}
			    }
			}
		}
	}
    }

    public Parameters createParameters(String line) {
	System.out.println("xxxxx "+line);
	return new DistanceConstrainParameters(line);
    }

    public Parameters parameters(Object obj) {
	throw new RuntimeException("This method has not been implemented yet");
    }


    private static class SatisfiedConstrainsFilter implements Filter {
	Protein model;
	double cutoff;
	public SatisfiedConstrainsFilter(Protein model, double cutoff) {
	    this.model = model;
	    this.cutoff = cutoff;
	}
	public String toString() {
	    return "filtering with "+cutoff;
	}
	
	public boolean accept(Object obj) {
	    DistanceConstrainParameters dcp = (DistanceConstrainParameters) obj;
	    Atom atom1 = model.getAtom(dcp.residue1Name(), dcp.residue1Number(), dcp.name1());
	    if (atom1 == null) throw new RuntimeException("cannot find "+dcp.residue1Name() +" "+
							 dcp.residue1Number()+" "+dcp.name1());
	    Atom atom2 = model.getAtom(dcp.residue2Name(), dcp.residue2Number(), dcp.name2()); 
	    if (atom2 == null) throw new RuntimeException("cannot find "+dcp.residue2Name()+" "+
							  dcp.residue2Number()+" "+dcp.name2());
	    if (Math.abs(2.0*(dcp.target-atom1.distanceFrom(atom2)))/(dcp.target+atom1.distanceFrom(atom2)) < cutoff) 
		return true;
	    else return false;
	}
    }
}
