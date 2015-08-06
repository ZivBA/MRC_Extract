package meshi.energy.simpleHydrogenBond;

import meshi.energy.AbstractEnergy;
import meshi.energy.TotalEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.hydrogenBond.AbstractHydrogenBond;
import meshi.geometry.hydrogenBond.AbstractHydrogenBondList;
import meshi.molecularElements.Atom;
import meshi.util.UpdateableException;

public class SimpleHydrogenBondEnergy extends AbstractEnergy {
	
	private AbstractHydrogenBondList hbList;
	private boolean toCalcDerivatives;
	private double hbEnergy;
	
	public SimpleHydrogenBondEnergy() {}

	public SimpleHydrogenBondEnergy(AbstractHydrogenBondList hbList, DistanceMatrix dm, double weight, boolean toCalcDerivatives) {
		super(toArray(dm, hbList) , weight);
		this.hbList = hbList; 
		this.toCalcDerivatives = toCalcDerivatives;
		comment = "HB";
	}
	
	
	/**
	 * The negative sum of the hbVal's from each hydrogen bond.
	 */
	public double evaluate() {
		if (!isOn())
			return 0;
		hbEnergy = 0.0;
		int hbListSize = hbList.size();
		AbstractHydrogenBond hb;
		for (int hbInd=0; hbInd<hbListSize ; hbInd++) {
			hb = hbList.getHB(hbInd);
			if (hb.active()) {
				hbEnergy -= hb.hbVal();
				if (toCalcDerivatives)
					hb.applyForcesToAtoms(-weight);
			}
		}
		return weight*hbEnergy;
	}

	/**
	 * The energy of each hydrogen bond is equally divided between its two polar heavy atoms. 
	 */
	public void evaluateAtoms() {
		double energy;
		int hbListSize = hbList.size();
		AbstractHydrogenBond hb;
		for (int hbInd=0; hbInd<hbListSize ; hbInd++) {
			hb = hbList.getHB(hbInd);
			if (hb.active()) {
				energy = -weight*hb.hbVal();
				hb.getFirstPolar().addEnergy(energy/2.0);
				hb.getSecondPolar().addEnergy(energy/2.0);
			}
		}
	}

	
	public void test(TotalEnergy totalEnergy, Atom atom) {
		if (hbList.size() == 0) 
			throw new RuntimeException("Cannot test "+this+"\n"+"No hydrogen bonds found");

		double[][] coordinates = new double[3][];
		coordinates[0] = atom.X();
		coordinates[1] = atom.Y();
		coordinates[2] = atom.Z();
		for(int i = 0; i< 3; i++) {
			try{totalEnergy.update();}catch(UpdateableException ue){}
			double x = coordinates[i][0];
			coordinates[i][1] = 0;
			double e1 = evaluate();
			double analiticalForce = coordinates[i][1];
			coordinates[i][0] += DX;
			// Whatever should be updated ( such as distance matrix torsion list etc. )
			try{totalEnergy.update();}catch(UpdateableException ue){}
			double e2 = evaluate();
			double de = e2-e1;
			double numericalForce = - de/DX;
			coordinates[i][0] -= DX;
			double diff = Math.abs(analiticalForce - numericalForce);

			if ((2*diff/(Math.abs(analiticalForce)+Math.abs(numericalForce)+verySmall)) > relativeDiffTolerance){
				System.out.println("Testing "+this);
				System.out.println("Atom["+atom.number()+"]."+XYZ.charAt(i)+" = "+x);
				System.out.println("Analytical force = "+analiticalForce);
				System.out.println("Numerical force  = "+numericalForce);

				System.out.println("diff = "+diff+"\n"+
						"tolerance = 2*diff/(|analiticalForce| + |numericalForce|+verySmall) = "+
						2*diff/(Math.abs(analiticalForce) + Math.abs(numericalForce)+verySmall));
				System.out.println();
			}
			if ((e1 == AbstractEnergy.INFINITY) | (e1 == AbstractEnergy.NaN))
				System.out.println("Testing "+this+"\ne1 = "+e1);
			if ((e2 == AbstractEnergy.INFINITY) | (e2 == AbstractEnergy.NaN))
				System.out.println("Testing "+this+"\ne2 = "+e2);
			if ((analiticalForce == AbstractEnergy.INFINITY) | (analiticalForce == AbstractEnergy.NaN))
				System.out.println("Testing "+this+"\nanaliticalForce = "+analiticalForce);

		}
	}
	
	public void printHBlist() {
		hbList.print();
	}
	
	public double hbEnergy() {
		return hbEnergy;
	}

}
