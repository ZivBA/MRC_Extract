package meshi.energy.hydrogenBondsPlane;

import meshi.geometry.Distance;
import meshi.molecularElements.Atom;

public class CNPlaneDistance {
private  final double K = 1;
private  double MINdis, MAXdis, slope, leftLineBound, rightLineBound;
private double a, b;
private double distanceVal, dDistanceDx, dDistanceDy, dDistanceDz;
private Atom  atom1, atom2;
private double energy;
protected double deDx1, deDy1, deDz1, deDx2, deDy2, deDz2;
private int swapFactor;

 protected CNPlaneDistance(double MINdis, double MAXdis, double slope) {
    this.MINdis = MINdis;
    this.MAXdis = MAXdis;
    this.slope = slope;
    leftLineBound = MINdis+slope;
    rightLineBound = MAXdis-slope;
 }

 protected  void set(Atom atom1, Atom atom2, Distance distance, int swapFactor) {
   this.atom1 = atom1;
   this.atom2 = atom2;
   this.swapFactor = swapFactor;
   distanceVal = distance.distance();
   dDistanceDx =  distance.dDistanceDx();
   dDistanceDy = distance.dDistanceDy();
   dDistanceDz = distance.dDistanceDz();
 }

    /**
    * energy and dirivarives calculation.
    **/
 protected  double  updateEnergy() {
    double d, d2, de;

    if (distanceVal <= MINdis || distanceVal >= MAXdis) {
            energy = 0;
            deDx1 = deDy1 = deDz1 = 0;
            deDx2 = deDy2 = deDz2 = 0;
            return energy;
    }
   if (distanceVal < leftLineBound)  {
       a = -2./(slope*slope*slope);
       b = 3./(slope*slope);
       d = distanceVal - MINdis;
       d2 = d*d;
       energy = K*(a*d2*d+b*d2);
       de = K*(3*a*d2+2*b*d);
   }
   else
       if (distanceVal > rightLineBound) {
       a = 2./(slope*slope*slope);
       b = -3./(slope*slope);
       d = distanceVal-rightLineBound;
       d2 = d*d;
       energy = K*(a*d2*d+b*d2+1);
       de = K*(3*a*d2+2*b*d);
       }
       else {
          energy = K;
          de = 0;
       }
   de = de*swapFactor;
   deDx1 = de*dDistanceDx;
   deDy1 = de*dDistanceDy;
   deDz1 = de*dDistanceDz;
   deDx2 = -deDx1;
   deDy2 = -deDy1;
   deDz2 = -deDz1;
 return energy;
 }

 protected Atom A() {return atom1;}
 protected Atom  B() {return atom2;}
 protected double distance(){return distanceVal;}
 protected double energy(){return energy;}

 public String toString() {
       return "\nDistanceElement:\n"+atom1+ "\n"+atom2+"\nenergy = "+energy+
       "\ndistance = "+distanceVal;
   }
}
