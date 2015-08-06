package meshi.energy.hydrogenBondsPlane;

import meshi.energy.NonBondedEnergyElement;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
/*       N1  --------- C2
 *          \         /
 *          C1 ---- N2
 */
/*
 * This class is an energy element for an Hydrogen Bond. Four atoms are involved: N1,C1 and N2,C2 atoms
 */
public class HydrogenBondsPlaneEnergyElement extends NonBondedEnergyElement {

    //---------------------------------- data filds -----------------------------------

    private DistanceMatrix distanceMatrix;

    private double weight;
    private double energy;
    private double energyDisC1N2, energyDisN1C2, energyTorsion;
    protected Atom C1, N1, C2, N2;
    protected CNtwoDistances twoDis;

    protected double force = -5;

   //---------fields needed for DistanceTerm----------------------
    protected static double MINdis = -10, MAXdis, slope = 1;
    protected final double DEFAULT_MAX_DIS = 5.5;
    protected Distance distanceN1C2, distanceC1N2;
    protected int swapFactorN1C2, swapFactorC1N2;
    private CNPlaneDistance distanceTermC1N2, distanceTermN1C2;
    //---------fields needed for Angle Term--------------------------
    private double maxAngleOfDistortion;
    private double slopeAngle;

    private CNPlaneTorsion torsionTerm;


    public HydrogenBondsPlaneEnergyElement(DistanceMatrix distanceMatrix,
                                           double weight, double maxAngleOfPlainDistortion) {
    	if (distanceMatrix.rMax()< DEFAULT_MAX_DIS)
        MAXdis = distanceMatrix.rMax();
        else MAXdis = DEFAULT_MAX_DIS;
        this.weight = weight;
        this.maxAngleOfDistortion = maxAngleOfPlainDistortion;
        slopeAngle = 1.57-maxAngleOfPlainDistortion;
        if (slopeAngle < 0)
            throw new RuntimeException("MaxAngleOfDistortion for HydrogenBondsPlane should be less than 90");
        this.distanceMatrix = distanceMatrix;
        torsionTerm = new CNPlaneTorsion(distanceMatrix, this.maxAngleOfDistortion, slopeAngle);
        distanceTermC1N2 = new  CNPlaneDistance(MINdis, MAXdis, slope);
        distanceTermN1C2 = new  CNPlaneDistance(MINdis, MAXdis, slope);
    }

    /*
     * @parm obj should be Distance
     */
    public void set(Object obj){
        this.twoDis = (CNtwoDistances)obj;
        distanceN1C2 = twoDis.distance1();
        distanceC1N2 = twoDis.distance2();
     Atom a1, a2, a3, a4;
        a1 = distanceN1C2.atom1();
        a2 = distanceN1C2.atom2();
        a3 = distanceC1N2.atom1();
        a4 = distanceC1N2.atom2();
        swapFactorN1C2 = 1;
        swapFactorC1N2 = 1;
        CN_AtomAttribute atom1_attribute = (CN_AtomAttribute) a1.getAttribute(CN_AtomAttribute.key);
        if (atom1_attribute.isN) {
            N1 = a1;
            C2 = a2;
        }
        else {
        if (!atom1_attribute.isC) System.out.println("Check  CNenergyElement.set - not C and N in the thirst pair");
            N1 = a2;
            C2 = a1;
            swapFactorN1C2 = -1;
        }

       atom1_attribute = (CN_AtomAttribute) a3.getAttribute(CN_AtomAttribute.key);
        if (atom1_attribute.isC) {
            C1 = a3;
            N2 = a4;
        }
        else {
        if (!atom1_attribute.isN) System.out.println("Check  CNenergyElement.set - not C and N in the second pair ");
            C1 = a4;
            N2 = a3;
            swapFactorC1N2 = -1;
        }
//        if (( C1.residueNumber() !=  N1.residueNumber()) || ( C2.residueNumber() !=  N2.residueNumber()))
        //System.out.println("HydrogenBondPlaneEnergyElement: Different residues!");
    }

    public void setAtoms(){
        atoms = new AtomList(4);
        atoms.fastAdd(N1);
        atoms.fastAdd(C2);
        atoms.fastAdd(N2);
        atoms.fastAdd(C2);
    }

     /*
    * use to evaluate and update some attribute of Distance to be use later by HydrogenBondsPairs
    */
    public double evaluate() {
            if (torsionTerm.set(N1,C1,N2,C2)) {
                distanceTermC1N2.set(C1,N2,distanceC1N2,swapFactorC1N2);
                distanceTermN1C2.set(N1,C2,distanceN1C2,swapFactorN1C2);
            energy = updateEnergy();
            updateAtoms();
            if ((! (energy < 0) ) & (! ( energy == 0)) & (!(energy > 0)))
                        System.out.println("weird energy "+energy);
            return energy*weight*force;
      }
       else return 0;
    }

   /**
    * energy and dirivarives calculation.
    **/
   public double updateEnergy() {
       energyDisC1N2 = distanceTermC1N2.updateEnergy();
       energyDisN1C2 = distanceTermN1C2.updateEnergy();
       energyTorsion = torsionTerm.updateEnergy();
       return energyDisC1N2*energyDisN1C2*energyTorsion;
   }

    public void updateAtoms(){
        double koef = -force*weight;
       if (! N1.frozen()) {
            N1.addToFx(koef* energyDisC1N2*
                   (torsionTerm.deDxN1*energyDisN1C2+energyTorsion*distanceTermN1C2.deDx1));
            N1.addToFy(koef* energyDisC1N2*
                   (torsionTerm.deDyN1*energyDisN1C2+energyTorsion*distanceTermN1C2.deDy1));
            N1.addToFz(koef* energyDisC1N2*
                   (torsionTerm.deDzN1*energyDisN1C2+energyTorsion*distanceTermN1C2.deDz1));
       }
        if (! C1.frozen()) {
            C1.addToFx(koef* energyDisN1C2*
                    (torsionTerm.deDxC1*energyDisC1N2+energyTorsion*distanceTermC1N2.deDx1));
            C1.addToFy(koef* energyDisN1C2*
                    (torsionTerm.deDyC1*energyDisC1N2+energyTorsion*distanceTermC1N2.deDy1));
            C1.addToFz(koef* energyDisN1C2*
                    (torsionTerm.deDzC1*energyDisC1N2+energyTorsion*distanceTermC1N2.deDz1));
        }
        if (! N2.frozen()) {
            N2.addToFx(koef*energyDisN1C2*
                    (torsionTerm.deDxN2*energyDisC1N2+energyTorsion*distanceTermC1N2.deDx2));
            N2.addToFy(koef*energyDisN1C2*
                    (torsionTerm.deDyN2*energyDisC1N2+energyTorsion*distanceTermC1N2.deDy2));
            N2.addToFz(koef* energyDisN1C2*
                    (torsionTerm.deDzN2*energyDisC1N2+energyTorsion*distanceTermC1N2.deDz2));
        }
         if (! C2.frozen()) {
            C2.addToFx(koef* energyDisC1N2*
                   (torsionTerm.deDxC2*energyDisN1C2+energyTorsion*distanceTermN1C2.deDx2));
            C2.addToFy(koef*energyDisC1N2*
                   (torsionTerm.deDyC2*energyDisN1C2+energyTorsion*distanceTermN1C2.deDy2));
            C2.addToFz(koef* energyDisC1N2*
                   (torsionTerm.deDzC2*energyDisN1C2+energyTorsion*distanceTermN1C2.deDz2));
         }
   }

   public final double energy(){return energy;}
   public final Atom N1() {return N1;}
   public final Atom C1() {return C1;}
   public final Atom N2() {return N2;}
   public final Atom C2() {return C2;}

   public String toString() {
        if ((C1 == null) & (N1 == null) & (C2 == null) & (N2 == null)) return "HydrogenBondEnergyElement: - atoms not set yet";
        if ((C1 == null) || (N1 == null) || (C2 == null) || (N2 == null)) throw new RuntimeException("HydrogenBondPlaneEnergyElement: This is weird\n"+torsionTerm.toString());
        return "HydrogenBondsPlaneEnergyElement (C, N): \n"
                +distanceTermC1N2.toString()
                +distanceTermN1C2.toString()
                +torsionTerm.toString();
    }
}
