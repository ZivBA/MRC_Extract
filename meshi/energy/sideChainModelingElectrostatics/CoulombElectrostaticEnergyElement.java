package meshi.energy.sideChainModelingElectrostatics;

import meshi.energy.NonBondedEnergyElement;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;

/**
 *  This object represents Electrostatic Energy between two Atoms. 
 *  It updates the energy according to the change in their distance and position.
 **/

public  class CoulombElectrostaticEnergyElement extends NonBondedEnergyElement {
	
    public static final double MAX_ENERGY = 100;
    public static final double ALPHA = 0.5;
    protected DistanceMatrix distanceMatrix;
    protected Atom atom1, atom2;
    protected int atom1Number, atom2Number;
    protected double dielectricConstant, q1, q2 ;//q1 = charge of atom # 1, etc.
    protected boolean frozen;
    protected double energy; 
    protected double weight; 
    protected double rMax; // The maximum distance between two atoms
    //that will be considered when calculating the electrostatics energy. 
    protected double contact;
    private final int FIRST = 0,SECOND=1;
    protected ChargeParametersList parametersList; // Holds all the atoms' charges data.

    /**
     * default constructor
     *
     **/
    public  CoulombElectrostaticEnergyElement() {}
	
    /**
     * constructor
     * 
     * @param parametersList
     * @param distanceMatrix
     * @param weight
     * @param dielectricConstant
     **/
    public  CoulombElectrostaticEnergyElement(ChargeParametersList parametersList, DistanceMatrix distanceMatrix, 
                                              double weight, double dielectricConstant) {
        this.parametersList = parametersList;
        this.weight = weight;
        this.distanceMatrix = distanceMatrix;
        this.dielectricConstant = dielectricConstant;
        rMax = distanceMatrix.rMax();
    }
	
    /**
     * setAtoms
     **/
    protected void setAtoms(){
        throw new RuntimeException("setAtoms() may not be used by ElectrostaticEnergyElement for efficiency.");
    }
    
    // 
    /**
     * Sets the relevant charge values for each atom in an atom pair.
     * The data is taken from the parameter list.
     * @param obj  an AtomPair
     **/
    public void set(Object obj) {
        Distance nonBonded = (Distance)obj;
        atoms = nonBonded.atoms();
        atom1 = nonBonded.atom1();
        atom2 = nonBonded.atom2();
        atom1Number = atom1.number();
        atom2Number = atom2.number();
        ChargeParameter parameter_1 = (ChargeParameter) parametersList.parameters(atom1);
        ChargeParameter parameter_2 = (ChargeParameter) parametersList.parameters(atom2);
        q1 = parameter_1.charge();
        q2 = parameter_2.charge();
    }
   
    /**
     * evaluate -
     * 1) Invokes updateEnergy() - updates the energy,
     * 2) Invokes updateAtoms() - updates the atoms' position.
     * @return double - energy*weight
     **/	
    public double evaluate() {
        updateEnergy();
        return energy*weight;
    }
    
    /**
     * Updates the energy
     * @return double - energy*weight
     **/	
    public double updateEnergy(){
        double EL;	
        double rMaxMinusDis;

        // if one of the atoms has zero charge, then the electostatic energy is 0.
        if(q1==0 || q2==0){
            energy = 0;
        }

        else {
            //double EL;	//electrostatics
            //double dELdD;
            double invD = 0;
            double multipleCharges = -99999; 
            double dis = -1;
			
            Distance distance;
			
            distance = distanceMatrix.distance(atom1Number, atom2Number);
            dis = distance.distance();
            rMaxMinusDis = rMax - dis;
            if (rMaxMinusDis <= 0) {
                energy = contact =  0;
                EL =0;
            }
            else {
            	invD = distance.invDistance();		
                multipleCharges = q1 * q2;
                EL =  multipleCharges * invD / dielectricConstant;
            	
/*
 * This is the original code in the current MESHI version. I changed it on 3/1/2008 so that
 * the dielectric constant will be a constant (1.0) and not 1/r  : 

            	invD = distance.invDistance();		
                invD2 = invD*invD;
                multipleCharges = q1 * q2;
                EL =  multipleCharges * invD2 / dielectricConstant;
*/		
            }



            //quench to zero in rMax
            double rMaxMinusDisSquare = rMaxMinusDis*rMaxMinusDis;
            double rMaxMinusDisSquarePlusAlpha = rMaxMinusDisSquare+ALPHA;
            contact = rMaxMinusDisSquare/rMaxMinusDisSquarePlusAlpha;
            energy = EL*contact;
				
        }
        return energy;
    }
   
   
    
    /**
     * toString
     * @return String
     **/
    public String toString() {
        if ((atom1 == null) & (atom2 == null)) return "ElectrostaticEnergyElement - atoms not set yet";
        if ((atom1 == null) | (atom2 == null)) throw new RuntimeException("This is weird\n"+
                                                                          "atom1 = "+atom1+"\n"+
                                                                          "atom2 = "+atom2); 
        Distance distance = distanceMatrix.distance(atom1Number, atom2Number);
        double dis = distance.distance();
        return ("ElectrostaticEnergyElement q1 = "+q1+" q2 = "+ q2 +" dielectricConstant = "+dielectricConstant+" Distance = "+
                dFormatSrt.f(dis)+" rMax = "+rMax+"\n"+atom1.verbose(1)+"\n"+atom2.verbose(1));
    }

}
