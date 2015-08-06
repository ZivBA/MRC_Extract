package meshi.energy.cAlphaPlane;

import meshi.geometry.Distance;
import meshi.molecularElements.Atom;

public class CAlphaPlaneDistance {
private  final double K = 1;
private  double MINdis, MAXdis, slope, leftLineBound, rightLineBound;   
private double a, b;
private double distanceAB, dDistanceDx, dDistanceDy, dDistanceDz;
private Atom  A, B; 
protected double energy;
protected double deDxA, deDyA, deDzA, deDxB, deDyB, deDzB;

 public CAlphaPlaneDistance(double MINdis, double MAXdis, double slope) {
    this.MINdis = MINdis;
    this.MAXdis = MAXdis;
    this.slope = slope;    
    leftLineBound = MINdis+slope;
    rightLineBound = MAXdis-slope;    
 }

 protected  void set(Atom A, Atom B, Distance distance ) {
   this.A = A;
   this.B = B;
   distanceAB = distance.distance();
   dDistanceDx = distance.dDistanceDx();
   dDistanceDy = distance.dDistanceDy();
   dDistanceDz = distance.dDistanceDz();   
 }   
 
    /**
    * energy and dirivarives calculation.
    **/
 protected  double  updateEnergy() {  
    double d, d2, de;
   if (distanceAB < leftLineBound)  {
       a = -2./(slope*slope*slope);
       b = 3./(slope*slope); 
       d = distanceAB - MINdis;
       d2 = d*d;
       energy = K*(a*d2*d+b*d2);       
       de = K*(3*a*d2+2*b*d);
   }
   else
       if (distanceAB > rightLineBound) {
       a = 2./(slope*slope*slope);
       b = -3./(slope*slope); 
       d = distanceAB-rightLineBound;
       d2 = d*d;
       energy = K*(a*d2*d+b*d2+1);       
       de = K*(3*a*d2+2*b*d);
       }
       else {
          energy = K;
          de = 0;
       }
   deDxA = de*dDistanceDx;
   deDyA = de*dDistanceDy;
   deDzA = de*dDistanceDz;
   deDxB = -deDxA;
   deDyB = -deDyA;
   deDzB = -deDzA;      
    return energy;
 }
  
 protected Atom A() {return A;}
 protected Atom  B() {return B;}   
 protected double distance(){return distanceAB;}
 protected double energy(){return energy;}
   
 public String toString() {
       return "DistanceElement:\n"+A+ "\n"+B+"\nenergyAB = "+energy+
       "\ndistanceAB = "+distanceAB;       
   }
 
}
