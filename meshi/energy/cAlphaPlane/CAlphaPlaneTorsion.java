package meshi.energy.cAlphaPlane;

import meshi.geometry.DisposableAngle;
import meshi.geometry.DisposableTorsion;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;



public class CAlphaPlaneTorsion {
/*       C              E
 *          \         /
 *           A ---- B
 *         /          \
 *    D                 F
 */
private final double MAXtorsion = 1; // |max(cos)|
private final double K = 1, a, c;
private final double angleNoDef; //parameters for Angle Term
private final double MINangle, MAXangle, slopeAng;//MINangle = angleNoDef MAXangle = Math.PI-angleNoDef
private DisposableAngle angleCAB, angleBAD, angleABE, angleABF;

private CAlphaPlaneAngle angleElementCAB, angleElementABE;
private Atom A, B, C, E;

private DisposableTorsion torsionCABE;
protected double dtorCABEeDxC, dtorCABEeDyC, dtorCABEeDzC, dtorCABEeDxA, dtorCABEeDyA, dtorCABEeDzA,
                 dtorCABEeDxB, dtorCABEeDyB, dtorCABEeDzB, dtorCABEeDxE, dtorCABEeDyE, dtorCABEeDzE;
protected double energyCABE;
protected double dCABEeDxC, dCABEeDyC, dCABEeDzC, dCABEeDxA, dCABEeDyA, dCABEeDzA,
                 dCABEeDxB, dCABEeDyB, dCABEeDzB, dCABEeDxE, dCABEeDyE, dCABEeDzE;

private DistanceMatrix distanceMatrix;
boolean smallAngleCABE;

public CAlphaPlaneTorsion(double angleNoDef, double slopeAng, DistanceMatrix distanceMatrix){
    this.slopeAng = slopeAng;
    this.angleNoDef = angleNoDef;
    MINangle = angleNoDef;
    MAXangle = Math.PI-angleNoDef;
    this.distanceMatrix = distanceMatrix;
    double MAXtorsion2 = MAXtorsion*MAXtorsion;
    a = -K/(MAXtorsion2*MAXtorsion2);
    c = 2*K/MAXtorsion2;
    angleElementCAB = new CAlphaPlaneAngle(MINangle, MAXangle, slopeAng);
    angleElementABE = new CAlphaPlaneAngle(MINangle, MAXangle, slopeAng);
}

  //to set and to check the space disposition conditions
  protected  void setTorsions(Atom A, Atom B, Atom C, Atom E){
  double angle;
    this.A = A;
    this.B = B;
    this.C = C;
    this.E = E;
    smallAngleCABE = false;

  //vectors - DC , AB
   angleCAB = new DisposableAngle(C,A,B,distanceMatrix);
   angle = angleCAB.angle();
   if ((Math.abs(angle) > MAXangle) || (Math.abs(angle) < MINangle)){
       smallAngleCABE = true;
   }

  //vectors - FE , BA
   angleABE = new DisposableAngle(A,B,E,distanceMatrix);
   angle = angleABE.angle();
   if ((Math.abs(angle) > MAXangle) || (Math.abs(angle) < MINangle)){
       smallAngleCABE = true;
   }
}

   /**
    * energy and dirivarives calculation.
    **/
 protected double updateEnergy() {
//angleElementBAC includes atoms B (atom1), A (atom2), C (atom3)
//angleElementABE includes atoms A, B, E
 if (smallAngleCABE) {
      energyCABE = 0;
      dCABEeDxA = dCABEeDyA = dCABEeDzA = 0;
      dCABEeDxB = dCABEeDyB = dCABEeDzB = 0;
      dCABEeDxC = dCABEeDyC = dCABEeDzC = 0;
      dCABEeDxE = dCABEeDyE = dCABEeDzE = 0;
  }
  else {
        	angleElementCAB.set(C,A,B,angleCAB);
            angleElementCAB.updateEnergy();
        	angleElementABE.set(A,B,E,angleABE);
        	angleElementABE.updateEnergy();

   	getTorsionEnergy();
  //torsion CABE, angles BAC and ABE
  dCABEeDxA = energyCABE*angleElementABE.energy()*angleElementCAB.deDx2 + //A is 2nd atom in the angleElementBAC
              energyCABE*angleElementCAB.energy()*angleElementABE.deDx1 + //A is 1st atom in the angleElementABE
              angleElementCAB.energy()*angleElementABE.energy()*dtorCABEeDxA;
  dCABEeDyA = energyCABE*angleElementABE.energy()*angleElementCAB.deDy2 + //A is 2nd atom in the angleElementBAC
              energyCABE*angleElementCAB.energy()*angleElementABE.deDy1 + //A is 1st atom in the angleElementABE
              angleElementCAB.energy()*angleElementABE.energy()*dtorCABEeDyA;
  dCABEeDzA = energyCABE*angleElementABE.energy()*angleElementCAB.deDz2 + //A is 2nd atom in the angleElementBAC
              energyCABE*angleElementCAB.energy()*angleElementABE.deDz1 + //A is 1st atom in the angleElementABE
              angleElementCAB.energy()*angleElementABE.energy()*dtorCABEeDzA;

  dCABEeDxB = energyCABE*angleElementABE.energy()*angleElementCAB.deDx3 + //B is 3d atom in the angleElementBAC
              energyCABE*angleElementCAB.energy()*angleElementABE.deDx2 + //B is 2nd atom in the angleElementABE
              angleElementCAB.energy()*angleElementABE.energy()*dtorCABEeDxB;
  dCABEeDyB = energyCABE*angleElementABE.energy()*angleElementCAB.deDy3 + //B is 3d atom in the angleElementBAC
              energyCABE*angleElementCAB.energy()*angleElementABE.deDy2 + //B is 2nd atom in the angleElementABE
              angleElementCAB.energy()*angleElementABE.energy()*dtorCABEeDyB;
  dCABEeDzB = energyCABE*angleElementABE.energy()*angleElementCAB.deDz3 + //B is 3d atom in the angleElementBAC
              energyCABE*angleElementCAB.energy()*angleElementABE.deDz2 + //B is 2nd atom in the angleElementABE
              angleElementCAB.energy()*angleElementABE.energy()*dtorCABEeDzB;

  dCABEeDxC = energyCABE*angleElementABE.energy()*angleElementCAB.deDx1 + //C is 1st atom in angleElementBAC
              0 + //C is not in angleElementABE
              angleElementCAB.energy()*angleElementABE.energy()*dtorCABEeDxC;
  dCABEeDyC = energyCABE*angleElementABE.energy()*angleElementCAB.deDy1 +//C is 1st atom in the angleElementBAC
              0 + //C is not in angleElementABE
              angleElementCAB.energy()*angleElementABE.energy()*dtorCABEeDyC;
  dCABEeDzC = energyCABE*angleElementABE.energy()*angleElementCAB.deDz1+//C is 1st atom in the angleElementBAC
              0 + //C is not in angleElementABE
              angleElementCAB.energy()*angleElementABE.energy()*dtorCABEeDzC;

  dCABEeDxE = 0 + //E is not in angleElementBAC
              energyCABE*angleElementCAB.energy()*angleElementABE.deDx3 + //E is 3d atom in the angleElementABE
              angleElementCAB.energy()*angleElementABE.energy()*dtorCABEeDxE;
  dCABEeDyE = 0 + //E is not in angleElementBAC
              energyCABE*angleElementCAB.energy()*angleElementABE.deDy3 + //E is 3d atom in the angleElementABE
              angleElementCAB.energy()*angleElementABE.energy()*dtorCABEeDyE;
  dCABEeDzE = 0 + //E is not in angleElementBAC
              energyCABE*angleElementCAB.energy()*angleElementABE.deDz3 + //E is 3d atom in the angleElementABE
              angleElementCAB.energy()*angleElementABE.energy()*dtorCABEeDzE;
  }
    return energyCABE*angleElementCAB.energy()*angleElementABE.energy();
}

 protected void getTorsionEnergy(){
//CABE
   if (smallAngleCABE) {
        energyCABE = 0;
        dtorCABEeDxA = dtorCABEeDyA = dtorCABEeDzA = 0;
        dtorCABEeDxB = dtorCABEeDyB = dtorCABEeDzB = 0;
        dtorCABEeDxC = dtorCABEeDyC = dtorCABEeDzC = 0;
        dtorCABEeDxE = dtorCABEeDyE = dtorCABEeDzE = 0;
   }
   else {
   double de;
   double cosTorsion;
   double cosTorsion2; // cosTorsion^2
   Atom atom1, atom2, atom3, atom4;
        torsionCABE = new DisposableTorsion (C,A,B,E, distanceMatrix);
        cosTorsion = torsionCABE.cosTorsion();
        cosTorsion2 = cosTorsion*cosTorsion;
	energyCABE = a*cosTorsion2*cosTorsion2+c*cosTorsion2;
	de = 4*a*cosTorsion2*cosTorsion+2*c*cosTorsion;
/*
    atom1 = torsionCABE.atom1;
	atom2 = torsionCABE.atom2;
	atom3 = torsionCABE.atom3;
	atom4 = torsionCABE.atom4;
        if (atom1 != C)
             System.out.println("CABE - not C1");
        if (atom2 != A)
             System.out.println("CABE - not A2");
        if (atom3 != B)
             System.out.println("CABE - not B3");
        if (atom4 != E)
             System.out.println("CABE - not D4\n");
//*/
//atom1 == C
        dtorCABEeDxC = de*torsionCABE.dCosTorsionDx1();
        dtorCABEeDyC = de*torsionCABE.dCosTorsionDy1();
        dtorCABEeDzC = de*torsionCABE.dCosTorsionDz1();
//atom2 == A
        dtorCABEeDxA = de*torsionCABE.dCosTorsionDx2();
        dtorCABEeDyA = de*torsionCABE.dCosTorsionDy2();
        dtorCABEeDzA = de*torsionCABE.dCosTorsionDz2();
//atom3 == B
        dtorCABEeDxB = de*torsionCABE.dCosTorsionDx3();
        dtorCABEeDyB = de*torsionCABE.dCosTorsionDy3();
        dtorCABEeDzB = de*torsionCABE.dCosTorsionDz3();
//atom4 == E
        dtorCABEeDxE = de*torsionCABE.dCosTorsionDx4();
        dtorCABEeDyE = de*torsionCABE.dCosTorsionDy4();
        dtorCABEeDzE = de*torsionCABE.dCosTorsionDz4();
   }
 }



public String toString() {
       return "TorsionElement:\nA:"+A+ "\nB:"+B+ "\nC:"+C+"\nE:"+E+
       "\nenergyCABE = "+energyCABE+
       "\ntorsionCABE.cos = "+torsionCABE.cosTorsion()+
       "\ntorsionCABE.angle = "+torsionCABE.torsion()*180/Math.PI+
       "\nangleElBAC = "+angleElementCAB+
       "\nangleElABE = "+angleElementABE;
   }
}

