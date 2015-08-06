package meshi.energy.distanceConstrains;
import java.util.Iterator;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomPair;
import meshi.molecularElements.Protein;
import meshi.util.MeshiList;
/**
 * Linearly constrains predefined pairs of atoms to a target distance.
**/
public class DistanceConstrainsEnergy extends SimpleEnergyTerm{
    private Protein protein;

    public DistanceConstrainsEnergy() {}
    

    public DistanceConstrainsEnergy(Protein protein, 
				    DistanceConstrainParametersList distanceConstrainParametersList, 
				    double weight) {
	
	super(toArray( ), distanceConstrainParametersList, weight);
	comment = "DistanceConstrains";
	this.weight = weight;
	this.protein = protein;
	createElementsList(distanceConstrainParametersList);
    }
    
    public void createElementsList(MeshiList list) {
	DistanceConstrainParametersList parametersList = (DistanceConstrainParametersList) list;
	elementsList = new MeshiList();
	DistanceConstrainElement distanceConstrain;
        DistanceConstrainParameters parameters;
	Iterator parametersIterator = parametersList.iterator();
	AtomPair constrainedAtoms;
	while ((parameters = (DistanceConstrainParameters)parametersIterator.next()) != null) {
	    constrainedAtoms = getAtomPair(protein, parameters);
	    elementsList.add(createElement(constrainedAtoms,parameters));
	}
    }
 
    public  EnergyElement createElement(Object obj, Parameters parameters) {
	AtomPair constrainedAtoms = (AtomPair) obj;
	return new DistanceConstrainElement( constrainedAtoms.atom1(), constrainedAtoms.atom2(), 
				      (DistanceConstrainParameters) parameters,this.weight);
    }

    private static AtomPair getAtomPair(Protein protein, DistanceConstrainParameters parameters) {
	Atom atom1 = protein.getAtom(parameters.residue1Name, parameters.residue1Number, parameters.name1);
	if (atom1 == null) throw new RuntimeException("cannot find "+parameters.residue1Name+" "+
						      parameters.residue1Number+" "+parameters.name1);
	Atom atom2 = protein.getAtom(parameters.residue2Name, parameters.residue2Number, parameters.name2); 
	if (atom2 == null) throw new RuntimeException("cannot find "+parameters.residue2Name+" "+
						      parameters.residue2Number+" "+parameters.name2);
	return new AtomPair(atom1, atom2);
    }
}
	
	
