/**  Forcefield parameters for Helium taken from Muller-Plathe et al. 
 * J. Chem. Phys. 98(12)9895-9904 (1993)
 * The original values are in KJoules and nano-meters 
 * (0.0848 and 0.228 respectively) and were transformed to 
 * Kcalories and Angstroms (0.0848/4.18 = 0.020287081 & 2.28)
 *
 **/
package meshi.molecularElements.helium;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;


public class Helium extends Atom implements HeliumInterface{
    public static final String comment = 
	"Forcefield parameters for Helium taken "+
	"from Muller-Plathe et al. "+
	"J. Chem. Phys. 98(12)9895-9904 (1993)";
    /**
     * A Helium atom in a specified position
     **/
    public Helium(double x, double y, double z) {
	super(x, y, z,
	      "HE",//name
	      "HEL",
	      1,
	      HE); //type
    }
	
    public Helium(double radius) {
	super(0,0,0,radius,
	      "HE",//name
	      "HEL",
	      1,
	      HEL); //type
   }
    public Helium(double radius, double minDistanceSqr, AtomList atomList) {
	super(0,0,0,radius, minDistanceSqr, 
	      "HE",//name
	      "HEL",
	      1,
	      HEL, atomList); //type
    }    
}
