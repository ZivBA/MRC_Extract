package meshi.geometry;
import meshi.util.Updateable;
import meshi.util.UpdateableException;
/**
 *This class defines a torsion pair for the use of the various two torsion energies.
 **/
public class TorsionPair implements Updateable {
    private int numberOfUpdates = 0;
	private Torsion torsion1, torsion2;

	public TorsionPair(Torsion torsion1, Torsion torsion2) {
	    this.torsion1 = torsion1;
 	    this.torsion2 = torsion2;
	}
	
	public Torsion torsion1() {return torsion1;}
	
	public Torsion torsion2() {return torsion2;}
     

    public void update(int numberOfUpdates) throws UpdateableException{
	if (numberOfUpdates == this.numberOfUpdates+1) {
	   torsion1.update(numberOfUpdates);
	   torsion2.update(numberOfUpdates);
	    this.numberOfUpdates++;
	}
	else if (numberOfUpdates != this.numberOfUpdates) 
	    throw new RuntimeException("Something weird with TorsionPair.update(int numberOfUpdates)\n"+
				       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
    }

        
    public String toString (){ 
    return "Torsion 1:\n-------------\n" + torsion1 + "\nTorsion 2:\n-------------\n" + torsion2 + "\n******************************\n\n";
}

}

