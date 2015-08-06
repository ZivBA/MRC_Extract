package meshi.energy.cAlphaHydrogenBond;
import meshi.energy.EnergyElement;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;

/*       C        E
 *        \       \
 *         A ---- B
 *       /       /
 *    D        F
 */
public class CAlphaEnergyElement extends EnergyElement {
private  final double weight;
private Atom A, B, C, D, E, F; 
protected double energy, energyDis, energyGeometry;

protected static final double MINdisCA = 4.5, MAXdisCA = 7.0, slope = 1;
protected static final double Q_H = -0.6;
protected final double MAXabsV2 = 13.4*13.4;

private Distance distanceAB;
private CAlphaDistance distanceElement; 
private CAlphaGeometry geometryElement;    
//private DistanceMatrix distanceMatrix;
private  ResidueList residues;
int resNumber;
  
public CAlphaEnergyElement(double weight, DistanceMatrix distanceMatrix,  ResidueList residues) {
    distanceElement = new CAlphaDistance(MINdisCA, MAXdisCA, slope);
    geometryElement = new CAlphaGeometry(MAXabsV2,distanceMatrix);      
    if (distanceMatrix.rMax() <  MAXdisCA)
        throw new RuntimeException("The current value of distanceMatrix.rMax = "+distanceMatrix.rMax()+
                " is too small for cAlphaHydrogenBond." +
                "\nAdvisable value of distanceMatrix.rMax should be no smaller than "+MAXdisCA);
    //this.distanceMatrix = distanceMatrix;
    this.weight = weight;
    this. residues = residues;
    resNumber = residues.size();
}

    protected boolean set(Distance distance){
    double distanceVal = distance.distance();
    if ((distanceVal < MINdisCA) || (distanceVal >  MAXdisCA)) return false;
        this.distanceAB = distance;
        Integer resNumberA, resNumberB;
        Residue resC, resE, resD, resF;
            A = distance.atom1();
        if (A.residue().dummy())
            return false;
            B = distance.atom2();
        if (B.residue().dummy())
            return false;
            resNumberA = A.residueNumber();
        if (resNumberA >= resNumber-1 || resNumberA == 0) return false;
            resNumberB = B.residueNumber();
        if (resNumberB >= resNumber-1 || resNumberB == 0) return false;
               resC = residues. residueAt(resNumberA-1);
               if (resC.dummy()) return false;
            resD = residues. residueAt(resNumberA+1);
                if (resD.dummy()) return false;
            resE = residues. residueAt(resNumberB-1);
               if (resE.dummy()) return false;
           resF = residues. residueAt(resNumberB+1);
               if (resF.dummy()) return false;
                C = resC.ca();
                D = resD.ca();
                F = resF.ca();
                E = resE.ca();
         //   System.out.println(distanceAB);
        //         System.out.println(A+"\n"+B+"\n");
        return true;
    }


 public double evaluate() {
        distanceElement.set(A, B, distanceAB);
        energyDis = distanceElement.updateEnergy();
        geometryElement.setGeometry(A,B,C,D,E,F,distanceAB);
        energyGeometry = geometryElement.updateEnergy();
        updateAtoms();
        return energyDis*energyGeometry*Q_H*weight;
   } 
 
protected void updateAtoms(){
double koef = Q_H*weight;        
double deDx, deDy, deDz;
   if (! A.frozen()) {
         deDx = distanceElement.deDxA*energyGeometry;
         deDx+= (energyDis*(geometryElement.dABCDeDxA+geometryElement.dBAEFeDxA));
         deDy =distanceElement.deDyA*energyGeometry;
         deDy += (energyDis*(geometryElement.dABCDeDyA+geometryElement.dBAEFeDyA));
         deDz =distanceElement.deDzA*energyGeometry;
         deDz += (energyDis*(geometryElement.dABCDeDzA+geometryElement.dBAEFeDzA));
         A.addToFx(-deDx*koef); 
         A.addToFy(-deDy*koef); 
         A.addToFz(-deDz*koef);          
        }
   if (!B.frozen()) {         
         deDx = distanceElement.deDxB*energyGeometry;
         deDx+= (energyDis*(geometryElement.dABCDeDxB+geometryElement.dBAEFeDxB));
         deDy =distanceElement.deDyB*energyGeometry;
         deDy += (energyDis*(geometryElement.dABCDeDyB+geometryElement.dBAEFeDyB));
         deDz =distanceElement.deDzB*energyGeometry;
         deDz += (energyDis*(geometryElement.dABCDeDzB+geometryElement.dBAEFeDzB));           
         B.addToFx(-deDx*koef); 
         B.addToFy(-deDy*koef); 
         B.addToFz(-deDz*koef);          
   }
   if (!C.frozen()) {         
        deDx = energyDis*geometryElement.dABCDeDxC;
        deDy = energyDis*geometryElement.dABCDeDyC;
        deDz = energyDis*geometryElement.dABCDeDzC;
        C.addToFx(-deDx*koef); 
        C.addToFy(-deDy*koef); 
        C.addToFz(-deDz*koef);          
   }

   if (!D.frozen()) {         
        deDx = energyDis*geometryElement.dABCDeDxD;
        deDy = energyDis*geometryElement.dABCDeDyD;
        deDz = energyDis*geometryElement.dABCDeDzD;
        D.addToFx(-deDx*koef); 
        D.addToFy(-deDy*koef); 
        D.addToFz(-deDz*koef);          
   }

   if (!E.frozen()) {         
        deDx = energyDis*geometryElement.dBAEFeDxE;
        deDy = energyDis*geometryElement.dBAEFeDyE;
        deDz = energyDis*geometryElement.dBAEFeDzE;
        E.addToFx(-deDx*koef); 
        E.addToFy(-deDy*koef); 
        E.addToFz(-deDz*koef);          
   }

   if (!F.frozen()) {         
        deDx = energyDis*geometryElement.dBAEFeDxF;
        deDy = energyDis*geometryElement.dBAEFeDyF;
        deDz = energyDis*geometryElement.dBAEFeDzF;
        F.addToFx(-deDx*koef); 
        F.addToFy(-deDy*koef); 
        F.addToFz(-deDz*koef);          
   }
}
    
  protected void setAtoms() {   
   atoms = new AtomList(6);
   atoms.add(A);
   atoms.add(B); 
   atoms.add(C);
   atoms.add(D); 
   atoms.add(E);
   atoms.add(F);     
   }
 
 public String toString() {
      return distanceElement+"\n"+geometryElement;
   }
    
}
