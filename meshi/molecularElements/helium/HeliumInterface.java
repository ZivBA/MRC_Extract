package meshi.molecularElements.helium;
import meshi.molecularElements.Atom;
import meshi.molecularElements.Residue;
public interface HeliumInterface {
    public static final int HE = Atom.addType("HE","DONE"); 
    public static final int HEL = Residue.addName("HEL","H",-1); 
    public static final String LENNARD_JONES_PARAMETERS = "meshi/parameters/helium/lennardJones.dat";
}
    
