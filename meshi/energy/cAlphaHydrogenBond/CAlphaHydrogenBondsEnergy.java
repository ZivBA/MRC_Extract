/*       C       E
 *        \       \
 *         A ---- B               
 *       /       /
 *     D        F       
 */

package meshi.energy.cAlphaHydrogenBond;

import java.util.Iterator;

import meshi.energy.AbstractEnergy;
import meshi.energy.TotalEnergy;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.Protein;
import meshi.parameters.AtomTypes;

public class CAlphaHydrogenBondsEnergy extends AbstractEnergy  implements AtomTypes {        

private  CAlphaHydrogenBondsParametersList parametersList;       
private CAlphaEnergyElement energyElement;
private DistanceList nonBondedList;//, bondedList;

public CAlphaHydrogenBondsEnergy(){super();}
     
 public CAlphaHydrogenBondsEnergy ( Protein protein, DistanceMatrix distanceMatrix,
    							CAlphaHydrogenBondsParametersList parametersList,								
								double weight){    	
    super(toArray(distanceMatrix),parametersList,weight);
    comment = "CalphaHB";
    energyElement = new CAlphaEnergyElement(weight,distanceMatrix,protein.residues());
    this.parametersList = parametersList;
    nonBondedList = distanceMatrix.nonBondedList();
 //   bondedList = distanceMatrix.bondedList();
 }	
 
 public double evaluate() { 
 if (! on) return 0.0;
 double energy = 0;
     Iterator nonBondedIter = nonBondedList.iterator();
     Distance pair;
     while ((pair  = (Distance) nonBondedIter.next()) != null) {
         if (!accept(pair)) continue;
        if (energyElement.set(pair))
            energy += energyElement.evaluate();
    }
/*
     Iterator bondedIter = bondedList.iterator();
     while ((pair  = (Distance) bondedIter.next()) != null) {
         if (!accept(pair)) continue;
        if (energyElement.set(pair)) {
            energy += energyElement.evaluate();
            energyElement.updateAtoms();
        }
    }    */

 return energy;
}
  
   public void evaluateAtoms() {
       throw new RuntimeException ("evaluateAtoms isn't used inside cAphaHydrogenBondsEnergy");
   } 
 
 	public boolean accept(Distance dis) {
            String atom1SS = dis.atom1().residue().secondaryStructure();
    		String atom2SS = dis.atom2().residue().secondaryStructure();
    		if (atom1SS == "COIL" || atom2SS == "COIL")
    			return false;
    		if (atom1SS == "HELIX" || atom2SS == "HELIX")
    			return false;
       int residueNumber1 = dis.atom1().residueNumber();
       int residueNumber2 = dis.atom2().residueNumber();
       return ((Math.abs(residueNumber1-residueNumber2) > 3));
  	}

    /* (non-Javadoc)
      * @see meshi.energy.AbstractEnergy#test(boolean, boolean)
      */
    public void test(TotalEnergy totalEnergy, Atom atom) {
          if (! on) System.out.println(""+this +" is off");
            Iterator nonBondedIter = nonBondedList.iterator();
            Distance pair;
            while ((pair  = (Distance) nonBondedIter.next()) != null) {
               if (!accept(pair)) continue;
               if (energyElement.set(pair)) {
                energyElement.setAtoms();
                energyElement.test(totalEnergy,atom);
               }
           }
/*            Iterator bondedIter = bondedList.iterator();
            while ((pair  = (Distance) bondedIter.next()) != null) {
                if (!accept(pair)) continue;
                if (energyElement.set(pair)) {
                    energyElement.setAtoms();
                    energyElement.test(totalEnergy,atom);
            }
          }
  */  }
}
