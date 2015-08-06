package meshi.energy.torsionVal;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.geometry.Angle;
import meshi.geometry.Torsion;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;

public class TorsionValEnergyElement extends EnergyElement {

    protected Atom atom1, atom2, atom3, atom4;
    protected double target;
    protected Torsion torsion;
    protected double weight;
    private int torCode;
    private boolean targetSet = false;

    public TorsionValEnergyElement(Torsion torsion, Parameters parameters, double weight) {
    	this.torsion = torsion;
    	target = torsion.torsion();
    	if (!((torsion.getTorsionName().equals("CHI1"))||
    			(torsion.getTorsionName().equals("CHI2"))||
    			(torsion.getTorsionName().equals("CHI3"))||
    			(torsion.getTorsionName().equals("CHI4"))||
    			(torsion.getTorsionName().equals("PHI"))||
    			(torsion.getTorsionName().equals("PSI"))))
    		throw new RuntimeException("I can not handle torsions other than Phi, Psi or Chi1,2,3,4");
    	this.weight = weight;
    	atom1 = torsion.atom1;
    	atom2 = torsion.atom2;
    	atom3 = torsion.atom3;
    	atom4 = torsion.atom4;
    	setAtoms();

    	updateFrozen();
    }

    public Torsion torsion() {return torsion;}

    public int getTorCode() {return torCode;}

    public void setTarget(double target) {
    	this.target = target;
    	targetSet = true;
    }

    public void setAtoms() {
    	atoms = new AtomList();
    	atoms.add(atom1);
    	atoms.add(atom2);
    	atoms.add(atom3);
    	atoms.add(atom4);
    }

    public double evaluate() {
    	if (!targetSet) {
    		return 0.0;
    	}
    	double dEdTorsion;
    	double tmp,diff;
    	double energy = 0;

    	diff = torsion.torsion()-target;
    	if (diff<-Math.PI)
    		diff = 2*Math.PI + diff;
    	else if (diff>Math.PI)
    		diff = diff - 2*Math.PI;

    	tmp = 1 - Math.cos(diff);
    	energy = weight*tmp;
    	dEdTorsion = -weight*Math.sin(diff);


    	if (! atom1.frozen()) {
    		atom1.addToFx(dEdTorsion*torsion.dTorsionDx1());  
    		atom1.addToFy(dEdTorsion*torsion.dTorsionDy1());    
    		atom1.addToFz(dEdTorsion*torsion.dTorsionDz1());
    	}
    	if (! atom2.frozen()) {
    		atom2.addToFx(dEdTorsion*torsion.dTorsionDx2());   
    		atom2.addToFy(dEdTorsion*torsion.dTorsionDy2());    
    		atom2.addToFz(dEdTorsion*torsion.dTorsionDz2());    
    	}
    	if (! atom3.frozen()) {
    		atom3.addToFx(dEdTorsion*torsion.dTorsionDx3()); 
    		atom3.addToFy(dEdTorsion*torsion.dTorsionDy3());    
    		atom3.addToFz(dEdTorsion*torsion.dTorsionDz3());    
    	}
    	if (! atom4.frozen()) {
    		atom4.addToFx(dEdTorsion*torsion.dTorsionDx4()); 
    		atom4.addToFy(dEdTorsion*torsion.dTorsionDy4());    
    		atom4.addToFz(dEdTorsion*torsion.dTorsionDz4());    
    	}

    	return energy;
    }

    public String toString() {
    	return ("TorsionValEnergyElement "+torsion.name()+" target = "+Angle.rad2deg(target)+
    			" force = "+dFormatSrt.f(weight)+" torsion = "+
    			dFormatSrt.f(Angle.rad2deg(torsion.torsion()))+" energy = "+dFormatSrt.f(evaluate())+"\n"+
    			atom1.verbose(1)+"\n"+atom2.verbose(1)+"\n"+atom3.verbose(1)+"\n"+atom4.verbose(1));
    }
}
