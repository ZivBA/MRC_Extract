package meshi.energy.cAlphaHydrogenBond;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;

public class CAlphaGeometry{
/*       C        E
 *        \       \
 *         A ---- B
 *       /       /
 *    D        F
 */
private final double  MAXabsV2;
private final double K = 1, KdivMAXabsV4;
private Atom A, B, C, D, E, F;
protected double energyABCD, energyBAEF; 
protected double dABCDeDxA, dABCDeDyA, dABCDeDzA, dABCDeDxB, dABCDeDyB, dABCDeDzB,
                 dABCDeDxC, dABCDeDyC, dABCDeDzC, dABCDeDxD, dABCDeDyD, dABCDeDzD;
protected double dBAEFeDxA, dBAEFeDyA, dBAEFeDzA, dBAEFeDxB, dBAEFeDyB, dBAEFeDzB,
                 dBAEFeDxE, dBAEFeDyE, dBAEFeDzE, dBAEFeDxF, dBAEFeDyF, dBAEFeDzF;
  
private double multipABCD, multipBAEF;
private Distance distanceAB, distanceDC, distanceFE; 
private double vAB_x, vAB_y, vAB_z;
private double vDC_x, vDC_y, vDC_z;
private double vBA_x, vBA_y, vBA_z;
private double vFE_x, vFE_y, vFE_z;

private DistanceMatrix distanceMatrix;
private final double INFINITE_DISTANCE = Distance.INFINITE_DISTANCE; 
   
public CAlphaGeometry(double  MAXabsV2, DistanceMatrix distanceMatrix)
 {                  
 this.MAXabsV2 = MAXabsV2;     
 KdivMAXabsV4 = K/(MAXabsV2*MAXabsV2);
 this.distanceMatrix = distanceMatrix;
}       
  
   //checking of  space disposition conditions     
  protected  void setGeometry(Atom A, Atom B, Atom C, Atom D, Atom E, Atom F,Distance distanceAB){
    this.A = A;
    this.B = B;
    this.C = C;
    this.D = D;   
    this.E = E;
    this.F = F;      
  this.distanceAB = distanceAB;
   vAB_x = distanceAB.dx();
   vAB_y = distanceAB.dy();
   vAB_z = distanceAB.dz();        
//vectors - DC , AB
  distanceDC = distanceMatrix.distance(D,C);
  if (distanceDC.distance() == INFINITE_DISTANCE)
       distanceDC = new Distance(D,C);
   vDC_x = distanceDC.dx();
   vDC_y = distanceDC.dy();
   vDC_z = distanceDC.dz();
  multipABCD = vAB_x*vDC_x+vAB_y*vDC_y+vAB_z*vDC_z;
  // 2 vecors FE and BA   
  vBA_x = -vAB_x;
  vBA_y = -vAB_y;
  vBA_z = -vAB_z;
  distanceFE = distanceMatrix.distance(F,E);
  if (distanceFE.distance() == INFINITE_DISTANCE)
      distanceFE = new Distance(F,E);
   vFE_x = distanceFE.dx();
   vFE_y = distanceFE.dy();
   vFE_z = distanceFE.dz();      
      
   multipBAEF = vBA_x*vFE_x+vBA_y*vFE_y+vBA_z*vFE_z;      
 }
  
   /**
    * energy and dirivarives calculation.
    **/     
protected double updateEnergy() {    
  double de, x2, delta;
  double dxA, dyA, dzA, dxB, dyB, dzB, dxC, dyC, dzC, dxD, dyD, dzD, dxE, dyE, dzE, dxF, dyF, dzF;
   dxA =  vDC_x;  
   dyA =  vDC_y;
   dzA =  vDC_z;
   dxB = -vDC_x;
   dyB = -vDC_y;
   dzB = -vDC_z;
   dxC = -vAB_x;
   dyC = -vAB_y;
   dzC = -vAB_z;
   dxD =  vAB_x;
   dyD =  vAB_y;
   dzD =  vAB_z;        

    x2 = multipABCD*multipABCD;
  if (x2 >  MAXabsV2) {
       energyABCD = 0;
       de = 0;
  }
  else {             
            delta = x2 - MAXabsV2;                               
	    energyABCD = KdivMAXabsV4*delta*delta;
	    de = 4*KdivMAXabsV4*delta*multipABCD;
  }
   dABCDeDxA = de*dxA;
   dABCDeDyA = de*dyA;
   dABCDeDzA = de*dzA;
   dABCDeDxB = de*dxB;
   dABCDeDyB = de*dyB;
   dABCDeDzB = de*dzB;   
   dABCDeDxC = de*dxC;
   dABCDeDyC = de*dyC;
   dABCDeDzC = de*dzC;
   dABCDeDxD = de*dxD;
   dABCDeDyD = de*dyD;
   dABCDeDzD = de*dzD;      

//BAEF   
   dxA = -vFE_x;
   dyA = -vFE_y;
   dzA = -vFE_z;
   dxB =  vFE_x;   
   dyB =  vFE_y;   
   dzB =  vFE_z;      
   dxE = -vBA_x;
   dyE = -vBA_y;
   dzE = -vBA_z;
   dxF =  vBA_x;
   dyF =  vBA_y;
   dzF =  vBA_z;      
   x2 = multipBAEF*multipBAEF;
if (x2 > MAXabsV2) {
       energyBAEF = 0;
       de = 0;
 } 
      else {  
            delta = x2 - MAXabsV2;                               
	    energyBAEF = KdivMAXabsV4*delta*delta;
	    de = 4*KdivMAXabsV4*delta*multipBAEF;
 }
   dBAEFeDxA = de*dxA;
   dBAEFeDyA = de*dyA;
   dBAEFeDzA = de*dzA;
   dBAEFeDxB = de*dxB;
   dBAEFeDyB = de*dyB;
   dBAEFeDzB = de*dzB;   
   dBAEFeDxE = de*dxE;
   dBAEFeDyE = de*dyE;
   dBAEFeDzE = de*dzE;
   dBAEFeDxF = de*dxF;
   dBAEFeDyF = de*dyF;
   dBAEFeDzF = de*dzF;        
   return energyABCD + energyBAEF;
 }

protected double energyABCD(){return energyABCD;} //*koef
protected double energyBAEF(){return energyBAEF;} //*koef

public String toString() {
       return "CAlphaGeometryElement:\nA:"+A+ "\nB:"+B+ "\nC:"+C+ "\nD:"+D+ "\nE:"+E+ "\nF:"+F+"\nenergyABCD = "+energyABCD+"  energyBAEF = "+energyBAEF+
       "\ndistanceAB = "+distanceAB + "\ndistanceDC = "+distanceDC+ "\ndistanceFE = "+distanceFE+
       "\nmultABCD = "+multipABCD+"\nmultBAEF = "+multipBAEF;
   }

}

