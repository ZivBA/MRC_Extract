/**
 * An energy term that constrain all atoms to a sphere.
 **/
package meshi.energy.volumeConstrain;
import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.DummyResidue;
import meshi.util.SortableMeshiList;
public class VolumeConstrain extends SimpleEnergyTerm{
    protected double radius, energy;
    protected Atom center; 
    protected AtomList atomList;
    protected VolumeConstrainElement energyElement;
    public VolumeConstrain(AtomList atomList,
			   double centerx,  double centery,  double centerz,  
			   double radius, double weight) {
	super(toArray( ), new VolumeConstrainParametersList(),weight);
	this.atomList = atomList;
	this.radius = radius;
	center = new DummyAtom(centerx, centery, centerz);
	comment = "VolumeConstrain";
	createElementsList();
    }

    public void createElementsList() {
	createElementsList(atomList);
    } 

    public Parameters getParameters(Object baseElement, SortableMeshiList parametersList) {
	return new Parameters();
    }


    public boolean sortableParameters() { return false;}
    public Parameters createParameters(String s) {return null;}
    public Parameters getKey(Object o) {return null;}
    public EnergyElement createElement(Object baseElement, Parameters parameters){
	return new VolumeConstrainElement(center, (Atom) baseElement, radius, weight);
    }
    
    private static class DummyAtom extends Atom { 
	public static final String COMMENT = "DummyAtom";
	
	public DummyAtom(double centerx,  double centery,  double centerz) {
	    super(centerx,centery,centerz,
		  "DA",//name
		  (new DummyResidue(1)),// unspecified monomer object
		  -1);   // unspecified atom type
	}
	public String test() {
	    return "dummy atom";
	}
	public String comment() {
	    return COMMENT;
	}
    }

}
