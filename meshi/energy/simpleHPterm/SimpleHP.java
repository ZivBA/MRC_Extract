package meshi.energy.simpleHPterm;

import java.util.Iterator;

import meshi.energy.AbstractEnergy;
import meshi.energy.TotalEnergy;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.util.UpdateableException;
import meshi.util.mathTools.HB_Sigma;

public class SimpleHP extends AbstractEnergy {
	
	private final double TRANSITION_REGION = 0.25; // Angs
	private DistanceMatrix dm;
	private double weightHydrophobic;
	private double weightHydrophilic;
	private boolean toCalcDerivatives;
	private double hydrophobicEnergy;
	private double hydrophilicEnergy;
	private HB_Sigma hydrophilicSigma;
	private HB_Sigma hydrophobicSigma;
	
    public SimpleHP(){}

    public SimpleHP(DistanceMatrix distanceMatrix, 
    		double weightHydrophobic, 
    		double weightHydrophilic, 
    		double radiusHydrophobic, 
    		double radiusHydrophilic, 
    		boolean toCalcDerivatives) {
    	super(toArray(distanceMatrix), weightHydrophobic);
    	comment = "SimpleHP";
    	dm = distanceMatrix;
		this.weightHydrophobic = weightHydrophobic;
		this.weightHydrophilic = weightHydrophilic;
		this.toCalcDerivatives = toCalcDerivatives;
		hydrophilicSigma = new HB_Sigma(0.1,radiusHydrophilic-TRANSITION_REGION,radiusHydrophilic+TRANSITION_REGION,
				radiusHydrophilic+2*TRANSITION_REGION,0.95,0.05);
		hydrophobicSigma = new HB_Sigma(0.1,radiusHydrophobic-TRANSITION_REGION,radiusHydrophobic+TRANSITION_REGION,
				radiusHydrophobic+2*TRANSITION_REGION,0.95,0.05);
    }

	public double evaluate() {
		if (!isOn())
			return 0.0;
		hydrophilicEnergy = hydrophobicEnergy = 0.0;
		Iterator iter = dm.nonBondedList().iterator();
		Distance dis;
		Atom atom1, atom2;
		double tmpSigmaDeriv;
		while((dis = (Distance) iter.next()) != null) {
			atom1 = dis.atom1();
			atom2 = dis.atom2();
			if (!atom1.isBackbone && atom1.active())
				if (atom1.isCarbon) { // hydrophobic
					hydrophobicEnergy -= hydrophobicSigma.s(dis.distance());
					if (toCalcDerivatives) {
						tmpSigmaDeriv = -weightHydrophobic*hydrophobicSigma.s_tag(dis.distance());  
						if (!atom1.frozen()) {
							atom1.addToFx(-tmpSigmaDeriv*dis.dDistanceDx());    // force = -derivative
							atom1.addToFy(-tmpSigmaDeriv*dis.dDistanceDy());
							atom1.addToFz(-tmpSigmaDeriv*dis.dDistanceDz());
						}						
						if (!atom2.frozen()) {
							atom2.addToFx(tmpSigmaDeriv*dis.dDistanceDx());    
							atom2.addToFy(tmpSigmaDeriv*dis.dDistanceDy());
							atom2.addToFz(tmpSigmaDeriv*dis.dDistanceDz());
						}						
					}
				}
				else if (atom1.isNitrogen || atom1.isOxygen) { // hydrophilic
					hydrophilicEnergy += hydrophilicSigma.s(dis.distance());
					if (toCalcDerivatives) {
						tmpSigmaDeriv = weightHydrophilic*hydrophilicSigma.s_tag(dis.distance());   
						if (!atom1.frozen()) {
							atom1.addToFx(-tmpSigmaDeriv*dis.dDistanceDx());    // force = -derivative
							atom1.addToFy(-tmpSigmaDeriv*dis.dDistanceDy());
							atom1.addToFz(-tmpSigmaDeriv*dis.dDistanceDz());
						}						
						if (!atom2.frozen()) {
							atom2.addToFx(tmpSigmaDeriv*dis.dDistanceDx());    
							atom2.addToFy(tmpSigmaDeriv*dis.dDistanceDy());
							atom2.addToFz(tmpSigmaDeriv*dis.dDistanceDz());
						}						
					}					
				}
			if (!atom2.isBackbone && atom2.active())
				if (atom2.isCarbon) { // hydrophobic
					hydrophobicEnergy -= hydrophobicSigma.s(dis.distance());
					if (toCalcDerivatives) {
						tmpSigmaDeriv = -weightHydrophobic*hydrophobicSigma.s_tag(dis.distance());  
						if (!atom1.frozen()) {
							atom1.addToFx(-tmpSigmaDeriv*dis.dDistanceDx());    // force = -derivative
							atom1.addToFy(-tmpSigmaDeriv*dis.dDistanceDy());
							atom1.addToFz(-tmpSigmaDeriv*dis.dDistanceDz());
						}						
						if (!atom2.frozen()) {
							atom2.addToFx(tmpSigmaDeriv*dis.dDistanceDx());    
							atom2.addToFy(tmpSigmaDeriv*dis.dDistanceDy());
							atom2.addToFz(tmpSigmaDeriv*dis.dDistanceDz());
						}						
					}
				}
				else if (atom2.isNitrogen || atom2.isOxygen) { // hydrophilic
					hydrophilicEnergy += hydrophilicSigma.s(dis.distance());
					if (toCalcDerivatives) {
						tmpSigmaDeriv = weightHydrophilic*hydrophilicSigma.s_tag(dis.distance());   
						if (!atom1.frozen()) {
							atom1.addToFx(-tmpSigmaDeriv*dis.dDistanceDx());    // force = -derivative
							atom1.addToFy(-tmpSigmaDeriv*dis.dDistanceDy());
							atom1.addToFz(-tmpSigmaDeriv*dis.dDistanceDz());
						}						
						if (!atom2.frozen()) {
							atom2.addToFx(tmpSigmaDeriv*dis.dDistanceDx());   
							atom2.addToFy(tmpSigmaDeriv*dis.dDistanceDy());
							atom2.addToFz(tmpSigmaDeriv*dis.dDistanceDz());
						}						
					}					
				}
		}	
		return (weightHydrophilic*hydrophilicEnergy + weightHydrophobic*hydrophobicEnergy);		
	}
	
	public double hydrophobicEnergy() {
		return hydrophobicEnergy;
	}

	public double hydrophilicEnergy() {
		return hydrophilicEnergy;
	}

	public void evaluateAtoms() {
		if (!isOn())
			return;
		Iterator iter = dm.nonBondedList().iterator();
		Distance dis;
		Atom atom1, atom2;
		while((dis = (Distance) iter.next()) != null) {
			atom1 = dis.atom1();
			atom2 = dis.atom2();
			if (!atom1.isBackbone && atom1.active())
				if (atom1.isCarbon) { // hydrophobic
					atom1.addEnergy(-weightHydrophobic*hydrophobicSigma.s(dis.distance()));
				}
				else if (atom1.isNitrogen || atom1.isOxygen) { // hydrophilic
					atom1.addEnergy(weightHydrophilic*hydrophilicSigma.s(dis.distance()));
				}
			if (!atom2.isBackbone && atom2.active())
				if (atom2.isCarbon) { // hydrophobic
					atom2.addEnergy(-weightHydrophobic*hydrophobicSigma.s(dis.distance()));
				}
				else if (atom2.isNitrogen || atom2.isOxygen) { // hydrophilic
					atom2.addEnergy(weightHydrophilic*hydrophilicSigma.s(dis.distance()));
				}
		}		
	}

	public void test(TotalEnergy totalEnergy, Atom atom) {
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

}