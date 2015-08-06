package meshi.energy.excludedVolumeImprovedDistance;
import meshi.energy.NonBondedEnergyElement;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;


public  class EVenergyElement extends NonBondedEnergyElement {
	protected DistanceMatrix distanceMatrix;
	protected Atom atom1, atom2;
	protected Distance distance;
	int type;
	protected double sigma,mini,A,B,C,D;
	protected boolean frozen;
	protected double dEdD;
	protected double dEdX;
	protected double dEdY;
	protected double dEdZ;
	protected double energy;
	protected double weight;
	protected double rMax;
	protected EVenergyParametersList parametersList;
	protected int fragStart=100000000,fragEnd=1000000000;

	public  EVenergyElement() {}
	public  EVenergyElement(EVenergyParametersList parametersList, DistanceMatrix distanceMatrix,
			int type, double weight) {
		this.parametersList = parametersList;
		this.type = type;
		this.weight = weight;
		this.distanceMatrix = distanceMatrix;
		this.rMax = distanceMatrix.rMax();
		//----------------------
	}
	protected void setAtoms(){
		throw new RuntimeException("setAtoms() may not be used by Excluded Volume for "+
		"efficiency.");
	}

	public void setFragLimits(int start , int end) {
		fragStart = start;
		fragEnd = end;
	}


	public void set(Object obj) {
		distance = (Distance) obj;
		atoms = distance.atoms();
		atom1 = distance.atom1();
		atom2 = distance.atom2();
		EVenergyParameters parameters = (EVenergyParameters) parametersList.parameters(distance);
		sigma = parameters.sigma;
		mini = parameters.mini;
		A = parameters.A;
		B = parameters.B;
		C = parameters.C;
		D = parameters.D;
		if (sigma>rMax) 
			throw new RuntimeException("SoftExcluded Volume: sigma="+sigma+" and it is larger " +
					"than rMax="+rMax);
	}

	public double evaluate() {
		double dEdD;
		double dis = distance.distance();
		if (type==0) {    // This functional form of the repulsive part is similar to that of Summa and Levitt (2007)
			if (dis>mini || (atom1.residueNumber()==atom2.residueNumber())) {
				energy = 0;
			}
			else {
				energy = weight*(A*dis*dis*dis + B*dis*dis + C*dis + D);
				dEdD = weight*(3*A*dis*dis + 2*B*dis + C);
/* The original non-derivable form. Maybe it will be useful in the future:
 				if (dis>sigma) {
					D3 = auxA*(dis-mini)*weight;
					energy = D3*(dis-mini);
					dEdD = 2*D3;
				}
				else if (dis>(sigma-0.5)) {
					energy = weight*(1.0 + 10.0*(sigma-dis));
					dEdD = -weight*10.0;
				}
				else {
					energy = weight*(-24.0 + 60.0*(sigma-dis));   //   = 6.0 + 60.0*(sigma-0.5-dis);
					dEdD = -weight*60.0;
				}
*/
				if (! atom1.frozen()) {
					atom1.addToFx(-1*dEdD*distance.dDistanceDx()); // force = -derivative   
					atom1.addToFy(-1*dEdD*distance.dDistanceDy());
					atom1.addToFz(-1*dEdD*distance.dDistanceDz());
				}
				if (! atom2.frozen()) {
					atom2.addToFx(dEdD*distance.dDistanceDx());
					atom2.addToFy(dEdD*distance.dDistanceDy());
					atom2.addToFz(dEdD*distance.dDistanceDz());
				}
			}
		}
		else
			throw new RuntimeException("An unknown excluded volume type");
		return energy;
	}


	public String toString() {
		if ((atom1 == null) & (atom2 == null)) return "ExcludedVolumeEnergyElement - atoms not set yet";
		if ((atom1 == null) | (atom2 == null)) throw new RuntimeException("This is weird\n"+
				"atom1 = "+atom1+"\n"+
				"atom2 = "+atom2); 
		double dis = distance.distance();
		return ("ExcludedVolumeEnergyElement sigma = "+sigma+" Distance = "+
				dFormatSrt.f(dis)+" rMax = "+rMax+"\n"+atom1.verbose(1)+"\n"+atom2.verbose(1));
	}
}
