package meshi.energy.hydrogenBondsPlane;

import meshi.geometry.DisposableAngle;
import meshi.geometry.DisposableTorsion;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;


public class CNPlaneTorsion {
/*       N1  --------- C2
 *          \         /
 *          C1 ---- N2
 */

        private Atom C1, N1, C2, N2;
        protected double deDxC1, deDyC1, deDzC1, deDxN1, deDyN1, deDzN1,
                        deDxC2, deDyC2, deDzC2, deDxN2, deDyN2, deDzN2;

        //---------fields needed for Angle Term--------------------------
        private final double undefAngle, slopeAngle;
        private final double MINangle, MAXangle;//MINangle = angleNoDef; MAXangle = Math.PI-angleNoDef
        private DisposableAngle angleN1C1N2, angleC1N2C2;
        private CNPlaneAngle angleElementN1C1N2, angleElementC1N2C2;

        //--------- fields needed for Torsion energy calculationm -------
        private final double MAXtorsion = 1; // |max(cos)|
        private final double K = 1,  a,  c;
        private DisposableTorsion torsionN1C1N2C2;
        private double energyTorsion, energy;
        private double dtorDxC1, dtorDyC1, dtorDzC1, dtorDxN1, dtorDyN1, dtorDzN1,
                         dtorDxC2, dtorDyC2, dtorDzC2, dtorDxN2, dtorDyN2, dtorDzN2;
        private DistanceMatrix distanceMatrix;

        protected CNPlaneTorsion(DistanceMatrix distanceMatrix, double undefAngle, double slopeAngle) {
            this.distanceMatrix = distanceMatrix;
            this.undefAngle = undefAngle;
            this.slopeAngle = slopeAngle;
            MINangle = undefAngle;
            MAXangle = Math.PI - undefAngle;
            double MAXtorsion2 = MAXtorsion*MAXtorsion;
            a = -K/(MAXtorsion2*MAXtorsion2);
            c = 2*K/MAXtorsion2;
            angleElementN1C1N2 = new CNPlaneAngle(MINangle, MAXangle, slopeAngle);
            angleElementC1N2C2 = new CNPlaneAngle(MINangle, MAXangle, slopeAngle);
        }


    protected  boolean set(Atom N1, Atom C1, Atom N2, Atom C2){
        this.C1 = C1;
        this.N1 = N1;
        this.C2 = C2;
        this.N2 = N2;
        double angle;
        angleN1C1N2 = new DisposableAngle (N1,C1,N2,distanceMatrix);
        angle = angleN1C1N2.angle();
        if ((Math.abs(angle) > MAXangle) || (Math.abs(angle) < MINangle))
              return false;
        angleC1N2C2 = new DisposableAngle (C1,N2,C2,distanceMatrix);
        angle = angleC1N2C2.angle();
        if ((Math.abs(angle) > MAXangle) || (Math.abs(angle) < MINangle))
            return false;
     return true;
     }

    protected double updateEnergy() {
        angleElementC1N2C2.set(C1,N2,C2,angleC1N2C2);
        angleElementN1C1N2.set(N1,C1,N2,angleN1C1N2);
        angleElementC1N2C2.updateEnergy();
        angleElementN1C1N2.updateEnergy();

        energyTorsion = torsionEnergy();
//torsion N1C1N2C2, angles N1C1N2 and C1N2C2
 deDxC1 = energyTorsion*angleElementN1C1N2.energy()*angleElementC1N2C2.deDx1 + //C1 is 1st atom in the angleElementC1N2C2
          energyTorsion*angleElementC1N2C2.energy()*angleElementN1C1N2.deDx2 + //C1 is 2nd atom in the angleElementN1C1N2
          angleElementC1N2C2.energy()*angleElementN1C1N2.energy()*dtorDxC1;
 deDyC1 = energyTorsion*angleElementN1C1N2.energy()*angleElementC1N2C2.deDy1 + //C1 is 1st atom in the angleElementC1N2C2
          energyTorsion*angleElementC1N2C2.energy()*angleElementN1C1N2.deDy2 + //C1 is 2nd atom in the angleElementN1C1N2
          angleElementC1N2C2.energy()*angleElementN1C1N2.energy()*dtorDyC1;
 deDzC1= energyTorsion*angleElementN1C1N2.energy()*angleElementC1N2C2.deDz1 + //C1 is 1st atom in the angleElementC1N2C2
          energyTorsion*angleElementC1N2C2.energy()*angleElementN1C1N2.deDz2 + //C1 is 2nd atom in the angleElementN1C1N2
          angleElementC1N2C2.energy()*angleElementN1C1N2.energy()*dtorDzC1;

 deDxN2 = energyTorsion*angleElementN1C1N2.energy()*angleElementC1N2C2.deDx2 + //N2 is 2nd atom in the angleElementC1N2C2
          energyTorsion*angleElementC1N2C2.energy()*angleElementN1C1N2.deDx3 + //N2 is 3d atom in the angleElementN1C1N2
          angleElementC1N2C2.energy()*angleElementN1C1N2.energy()*dtorDxN2;
 deDyN2 = energyTorsion*angleElementN1C1N2.energy()*angleElementC1N2C2.deDy2 + //N2 is 2nd atom in the angleElementC1N2C2
          energyTorsion*angleElementC1N2C2.energy()*angleElementN1C1N2.deDy3 + //N2 is 3d atom in the angleElementN1C1N2
          angleElementC1N2C2.energy()*angleElementN1C1N2.energy()*dtorDyN2;
 deDzN2 = energyTorsion*angleElementN1C1N2.energy()*angleElementC1N2C2.deDz2 + //N2 is 2nd atom in the angleElementC1N2C2
          energyTorsion*angleElementC1N2C2.energy()*angleElementN1C1N2.deDz3 + //N2 is 3d atom in the angleElementN1C1N2
          angleElementC1N2C2.energy()*angleElementN1C1N2.energy()*dtorDzN2;

 deDxC2  = energyTorsion*angleElementN1C1N2.energy()*angleElementC1N2C2.deDx3 + //C2 is 3d atom in angleElementC1N2C2
          0 + //C2 is not in angleElementN1C1N2
          angleElementC1N2C2.energy()*angleElementN1C1N2.energy()*dtorDxC2;
 deDyC2  = energyTorsion*angleElementN1C1N2.energy()*angleElementC1N2C2.deDy3 +//C2 is 3d atom in the angleElementC1N2C2
          0 + //C2 is not in angleElementN1C1N2
          angleElementC1N2C2.energy()*angleElementN1C1N2.energy()*dtorDyC2;
 deDzC2 = energyTorsion*angleElementN1C1N2.energy()*angleElementC1N2C2.deDz3+//C2 is 3d atom in the angleElementC1N2C2
          0 + //C2 is not in angleElementN1C1N2
          angleElementC1N2C2.energy()*angleElementN1C1N2.energy()*dtorDzC2;

 deDxN1 = 0 + //N1 is not in angleElementC1N2C2
          energyTorsion*angleElementC1N2C2.energy()*angleElementN1C1N2.deDx1 + //N1 is 1st atom in the angleElementN1C1N2
          angleElementC1N2C2.energy()*angleElementN1C1N2.energy()*dtorDxN1;
 deDyN1 = 0 + //N1 is not in angleElementC1N2C2
          energyTorsion*angleElementC1N2C2.energy()*angleElementN1C1N2.deDy1 + //N1 is 1st atom in the angleElementN1C1N2
          angleElementC1N2C2.energy()*angleElementN1C1N2.energy()*dtorDyN1;
 deDzN1 = 0 + //N1 is not in angleElementC1N2C2
          energyTorsion*angleElementC1N2C2.energy()*angleElementN1C1N2.deDz1 + //N1 is 1st atom in the angleElementN1C1N2
          angleElementC1N2C2.energy()*angleElementN1C1N2.energy()*dtorDzN1;

        energy =  energyTorsion*angleElementN1C1N2.energy()*angleElementC1N2C2.energy();
        return energy;
    }

   private double torsionEnergy(){
       double de, energy;
       double cosTorsion;
       double cosTorsion2; // cosTorsion^2
            torsionN1C1N2C2 = new DisposableTorsion (N1,C1,N2,C2, distanceMatrix);
            cosTorsion = torsionN1C1N2C2.cosTorsion();
            cosTorsion2 = cosTorsion*cosTorsion;
            energy = a*cosTorsion2*cosTorsion2+c*cosTorsion2;
            de = 4*a*cosTorsion2*cosTorsion+2*c*cosTorsion;
/* //for checking only
       Atom atom1, atom2, atom3, atom4;
            atom1 = torsionN1C1N2C2.atom1;
            atom2 = torsionN1C1N2C2.atom2;
            atom3 = torsionN1C1N2C2.atom3;
            atom4 = torsionN1C1N2C2.atom4;
            if (atom1 != N1)
                 System.out.println("N1C1N2C2 - not N1");
            if (atom2 != C1)
                 System.out.println("N1C1N2C2 - not C1");
            if (atom3 != N2)
                 System.out.println("N1C1N2C2 - not N2");
            if (atom4 != C2)
                 System.out.println("N1C1N2C2 - not C2");
//*/
//atom1 == N1
        dtorDxN1 = de*torsionN1C1N2C2.dCosTorsionDx1();
        dtorDyN1 = de*torsionN1C1N2C2.dCosTorsionDy1();
        dtorDzN1 = de*torsionN1C1N2C2.dCosTorsionDz1();
//atom2 == C1
        dtorDxC1 = de*torsionN1C1N2C2.dCosTorsionDx2();
        dtorDyC1 = de*torsionN1C1N2C2.dCosTorsionDy2();
        dtorDzC1 = de*torsionN1C1N2C2.dCosTorsionDz2();
//atom3 == N2
        dtorDxN2 = de*torsionN1C1N2C2.dCosTorsionDx3();
        dtorDyN2 = de*torsionN1C1N2C2.dCosTorsionDy3();
        dtorDzN2 = de*torsionN1C1N2C2.dCosTorsionDz3();
//atom4 == C2
       dtorDxC2 = de*torsionN1C1N2C2.dCosTorsionDx4();
       dtorDyC2 = de*torsionN1C1N2C2.dCosTorsionDy4();
       dtorDzC2 = de*torsionN1C1N2C2.dCosTorsionDz4();
     return energy;
     }

    protected double energy(){return energy;}
    protected Atom N1() {return N1;}
    protected Atom C1() {return C1;}
    protected Atom N2() {return N2;}
    protected Atom C2() {return C2;}

    public String toString() {
       if ((C1 == null) & (N1 == null) & (C2 == null) & (N2 == null)) return "TorsionTerm: - atoms have not set yet";
       if ((C1 == null) || (N1 == null) || (C2 == null) || (N2 == null)) throw new RuntimeException("TorsionTerm: This is weird\n"+
               "Atom C1 = "+C1+"\n"+
               "Atom N1 = "+N1+"\n"+
               "Atom C2 = "+C2+"\n"+
               "Atom N2 = "+N2+"\n");
       return "TorsionTerm (C, N): \n"+
               "ResidueNumbers: Res1 = "+ C1.residueNumber()+
               "Res2 = "+C2.residueNumber()+
               "Atom C1 = "+C1+"\n"+
               "Atom N1 = "+N1+"\n"+
               "Atom C2 = "+C2+"\n"+
               "Atom N2 = "+N2+"\n"+
               angleElementN1C1N2.toString()+
               angleElementC1N2C2.toString();
   }
}

