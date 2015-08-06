package meshi.energy.torsionSpaceMinimization;

import meshi.geometry.ViewAt;
import meshi.molecularElements.Atom;

public class TempCalcForthAtom {
	
	// Auxilary arrays
	private double[] moveTo = new double[4];
    double[][] M = new double[4][4] ;
    double[][] invM = new double[4][4] ;
	private double[] xyza = new double[4];
	
	
	public double[] getCoors(double bond,double angle,double tor,Atom A1,Atom A2,Atom A3) {
		double[] results = new double[3];
		localGetCoordinates(bond, angle, tor, A1, A2, A3);
		results[0] = xyza[0];
		results[1] = xyza[1];
		results[2] = xyza[2];
		return results;
	}
	
	public void setForthAtom(double bond,double angle,double tor,Atom A1,Atom A2,Atom A3,Atom A4) {
		localGetCoordinates(bond, angle, tor, A1, A2, A3);
		A4.setXYZ(xyza[0], xyza[1], xyza[2]);
	}
	
	// Puting the coordinates in "xyza"
	private void localGetCoordinates(double bond, double angle, double tor, Atom A1, Atom A2, Atom A3) {
        /*xyza should be of type double[4]*/
        moveTo[0] = bond*Math.sin(angle)*Math.sin(tor);
        moveTo[1] =bond*Math.sin(angle)*Math.cos(tor);
        moveTo[2] =  bond*Math. cos(angle);
        moveTo[3] = 1;

        ViewAt.transformToOrigin(M,invM,A3, A2,A1);
        // moveTo*invM
        for (int i = 0; i < 4; i++) {
        	xyza[i] = 0;
        	for(int j = 0; j < 4; j++) {
                xyza[i] += moveTo[j] * invM[j][i];
            }
        }		
	}
	
}
