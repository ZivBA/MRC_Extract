package meshi.energy.cAlphaPlane;
import meshi.energy.EnergyElement;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
/*
*       C        E
 *        \          \
 *         A ---- B
 */
public class CAlphaPlaneEnergyElement extends EnergyElement {
private  final double weight;
private Atom A, B, C, E;
private Distance distance;
private double energy, energyDis, energyTorAng;

protected static final double MINdisCA = 4.5, MAXdisCA = 7.0, slope = 1; //parameters for Distance term
private final double angleNoDef = 6*Math.PI/180.0, slopeAng = 0.5; //parameters for Angle Term
protected static final double Q_H = -0.6; // common coefficient

private CAlphaPlaneDistance distanceElement; 
private CAlphaPlaneTorsion torsionElement;    
//private DistanceMatrix distanceMatrix;
private  ResidueList residues;
  
public CAlphaPlaneEnergyElement (double weight, DistanceMatrix distanceMatrix,  ResidueList residues) {
    distanceElement = new CAlphaPlaneDistance(MINdisCA, MAXdisCA, slope);    
    torsionElement = new CAlphaPlaneTorsion(angleNoDef, slopeAng, distanceMatrix);      
    if (distanceMatrix.rMax() <  MAXdisCA)
        throw new RuntimeException("The current value of distanceMatrix.rMax = "+distanceMatrix.rMax()+
                " is too small for cAlphaPlane." +
                "\nAdvisable value of distanceMatrix.rMax should be no smaller than "+MAXdisCA);
    this.weight = weight;
  //  this.distanceMatrix = distanceMatrix;
    this. residues = residues;
}

    protected boolean set(Distance distance){
    double distanceVal = distance.distance();
    if ((distanceVal < MINdisCA) || (distanceVal >  MAXdisCA)) return false;
        this.distance = distance;
        Integer resNumberA, resNumberB;
        Residue resC, resE;
            A = distance.atom1();
            B = distance.atom2();
            resNumberA = A.residueNumber();
            resNumberB = B.residueNumber();
               resC = residues. residueAt(resNumberA-1);
               if (resC.dummy()) return false;
               resE = residues. residueAt(resNumberB+1);
               if (resE.dummy()) return false;
                C = resC.ca();
                E = resE.ca();
        return true;
    }

 public double evaluate() {
        distanceElement.set(A, B, distance);
        energyDis = distanceElement.updateEnergy();               
        torsionElement.setTorsions(A,B,C,E);
        energyTorAng = torsionElement.updateEnergy();  
        updateAtoms();             
        return energyDis*energyTorAng*Q_H*weight;
   }
 

protected void updateAtoms(){
double koef = Q_H*weight;        
double deDx, deDy, deDz;
   if (! A.frozen()) {         
         deDx = distanceElement.deDxA*energyTorAng +
                energyDis*torsionElement.dCABEeDxA;
         deDy = distanceElement.deDyA*energyTorAng +
                energyDis*torsionElement.dCABEeDyA;
         deDz = distanceElement.deDzA*energyTorAng +
                energyDis*torsionElement.dCABEeDzA;
         A.addToFx(-deDx*koef); 
         A.addToFy(-deDy*koef); 
         A.addToFz(-deDz*koef);          
        }
   if (!B.frozen()) {         
         deDx = distanceElement.deDxB*energyTorAng +
                energyDis*torsionElement.dCABEeDxB;
         deDy = distanceElement.deDyB*energyTorAng +
                energyDis*torsionElement.dCABEeDyB;
         deDz = distanceElement.deDzB*energyTorAng +
                energyDis*torsionElement.dCABEeDzB;
         B.addToFx(-deDx*koef); 
         B.addToFy(-deDy*koef); 
         B.addToFz(-deDz*koef);          
   }
   if (!C.frozen()) {         
        deDx = energyDis*(torsionElement.dCABEeDxC);
        deDy = energyDis*(torsionElement.dCABEeDyC);
        deDz = energyDis*(torsionElement.dCABEeDzC);
        C.addToFx(-deDx*koef); 
        C.addToFy(-deDy*koef); 
        C.addToFz(-deDz*koef);          
   }
   if (!E.frozen()) {         
        deDx = energyDis*(torsionElement.dCABEeDxE);
        deDy = energyDis*(torsionElement.dCABEeDyE);
        deDz = energyDis*(torsionElement.dCABEeDzE);
        E.addToFx(-deDx*koef); 
        E.addToFy(-deDy*koef); 
        E.addToFz(-deDz*koef);          
   }
}
    
  protected void setAtoms() {   
   atoms = new AtomList(4);
   atoms.add(A);
   atoms.add(B); 
   atoms.add(C);
   atoms.add(E);
   }
 
 public String toString() {
      return distanceElement+"\n"+torsionElement;
   }    
}
