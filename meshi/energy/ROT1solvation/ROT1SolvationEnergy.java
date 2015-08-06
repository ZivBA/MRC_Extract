package meshi.energy.ROT1solvation;

import java.util.Iterator;

import meshi.energy.CooperativeEnergyTerm;
import meshi.energy.ROT1solvation.parameters.AbstractROT1Parameters;
import meshi.energy.ROT1solvation.parameters.WeightRepresentativePolars;
import meshi.energy.ROT1solvation.parameters.WeightsHydrophobicMarkerBaysian;
import meshi.energy.ROT1solvation.parameters.WeightsHydrophobicSA;
import meshi.energy.ROT1solvation.parameters.WeightsPolarSA;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.mathTools.Sigma;


public final class ROT1SolvationEnergy extends CooperativeEnergyTerm 
implements WeightsHydrophobicSA, WeightsPolarSA, WeightRepresentativePolars , WeightsHydrophobicMarkerBaysian {

     
    /** 
     * These are fields for temporary array results that are needed in the evaluation stage.
     * They are declared as fields so that time will not be waisted on creating new 
     * instances of the arrays.
     **/
     private double[] energyVals;     
     private double[] dEnergyVals;     
     private double[] AtomSumSigmC;
     private double[] forceX;
     private double[] forceY;       
     private double[] forceZ;
       
    /** 
     * These fields are for general use in the class
     **/       
     private int[] lut; 	// The look-up table (lut) converts the atom internal number (field of Atom) to its index in the atom list given to the constructor.      
     private int atomListSize;
     private AbstractROT1Parameters parameters; // The instance of the parameter list object.
     private int[] atomTypes; // The atom types of each atom in the protein.
     private ROT1SolvationCalculation solvCalc = null; // The encapsulation of the solvate calculation
     private boolean toCalcDerivatives = true; // If 'false' derivatives will not be calculated.
	 private  boolean[] matchingRes;

     // parameters for the sigmoid
     private double p1,p2,valAtp1,valAtp2,end,cutoff;
          
     // Fields for seeing the hydrophobic and polar effects
     private double energyHydrophobicWeightedSA, energyPolarWeightedSA, energyRepresentativePolar, energyHydrophobicBaysian; 
     
     
    public ROT1SolvationEnergy() {}
    
    public ROT1SolvationEnergy(AtomList atomList, 
                    DistanceMatrix dm,
				    AbstractROT1Parameters parameters,
				    boolean toCalcDerivatives,
				    boolean[] matchingRes,
				    double weight) {
	super(toArray(dm),atomList, dm, null, weight);
	int c;
	int maxAtomNum=-1;
	comment = "ROT1solv";
	atomListSize = atomList.size();
	this.parameters = parameters;
	this.matchingRes = matchingRes;
	if (parameters == null)
	    throw new RuntimeException("The parameters object for this Solvatation term is NULL");
	solvCalc = new ROT1SolvationCalculation(parameters);
	this.toCalcDerivatives = toCalcDerivatives;
	p1 = parameters.p1();
	p2 = parameters.p2();
	valAtp1 = parameters.valAtp1();
	valAtp2 = parameters.valAtp2();
	end = parameters.end();
	cutoff = 0.5*(p1+p2);
	

    // Creating the auxilary arrays
    energyVals = new double[atomListSize];
    dEnergyVals = new double[atomListSize];
    AtomSumSigmC = new double[atomListSize];
    forceX = new double[atomListSize];
    forceY = new double[atomListSize];
    forceZ = new double[atomListSize];
    
    
	// Creating the lookup table for the atom numbers.
	// The table converts the atom internal number (field of Atom) to its index  
	// in the atom list given to the constructor.      
	for (c=0; c<atomListSize ; c++) {
	    if (atomList.atomAt(c).number() > maxAtomNum)
	    	maxAtomNum = atomList.atomAt(c).number();
	}
	lut = new int[maxAtomNum+1];
	for (c=0; c<maxAtomNum ; c++) {
	    lut[c] = -1;
	}
	for (c=0; c<atomListSize ; c++) {
	    lut[atomList.atomAt(c).number()] = c;
	}
	
	// Setting up the atom type array, that can associate each atom in the protein with its atom type.
	atomTypes = new int[atomListSize];
	for (c=0; c<atomListSize ; c++) 
		atomTypes[c] = atomList.atomAt(c).type;
										
    } // of the constructor
    

/*
    private final double[] scSA = {20, //a
        	40, //c
            45, //d
            80, //e
           115, //f
           0, //g
           100, //h
            80, //i
           105, //k
            75, //l
            90, //m
            55, //n
            50, //p
            70, //q
           100, //r
            20, //s
            45, //t
            60, //v
           160, //w
           135}; //y 
    private final double bbSA = 55.0/75.0;  // Should be 95,  but this is without hydrogen bonding in a SS;
*/
    
    public void evaluateAtoms() {
		evaluate(true);
    }
    
    public double evaluate() {
    	return evaluate(false);
    }

    public double evaluate(boolean evaluateAtoms) {
	if (! on) return 0.0;
	double energy = 0;
	energyHydrophobicWeightedSA = 0; 
	energyPolarWeightedSA = 0;
	energyRepresentativePolar = 0;
	energyHydrophobicBaysian = 0;
	int cc;        
	DistanceList dislist = dm.nonBondedList();
	Iterator<Distance> iter;
	Distance dis;
	ROT1SolvationDistanceAttribute sigmaValues; 
	int ind1,ind2;
	Atom atom;


	//Reseting the auxilary arrays and variables
	if (toCalcDerivatives) {
		for (cc=0 ; cc<atomListSize ; cc++) {
			AtomSumSigmC[cc] = 0; 
			forceX[cc] = forceY[cc] = forceZ[cc] = 0.0;
		} 
	}
	else {
		for (cc=0 ; cc<atomListSize ; cc++) {
			AtomSumSigmC[cc] = 0; 		
		}
	}
	
	// First pass over the non-bonded list
	iter = dislist.iterator();
	if (toCalcDerivatives) {
		while((dis = iter.next()) != null) {
			if (parameters.filterDisForRelevance(dis)) {
				if (dis.getAttribute(ROT1SolvationDistanceAttribute.SOLVATE_ROT1_ATTRIBUTE) == null) {
					sigmaValues = new ROT1SolvationDistanceAttribute();
					dis.addAttribute(sigmaValues);
					if (dis.frozen)
						updateSigmVals(dis);
				}
				if (!dis.frozen) 
					updateSigmVals(dis);
				sigmaValues = (ROT1SolvationDistanceAttribute) dis.getAttribute(ROT1SolvationDistanceAttribute.SOLVATE_ROT1_ATTRIBUTE);
				ind1 = lut[dis.atom1().number()];
				ind2 = lut[dis.atom2().number()];
				AtomSumSigmC[ind1] += sigmaValues.sigmCa1;
				AtomSumSigmC[ind2] += sigmaValues.sigmCa2;
			}
		}
	}
	else {
		while((dis = iter.next()) != null) {
			if (parameters.filterDisForRelevance(dis) && (dis.distance()<cutoff)) {
// This is the code for the attunated distance attempt with CBs and CAs only:				
//				ind1 = lut[dis.atom1().number()];
//				ind2 = lut[dis.atom2().number()];
//				double factor = 36.0/Math.max(3.0, dis.distance())*1.0/Math.max(3.0, dis.distance());
//				if (dis.atom2().name().equals("CA"))
//					AtomSumSigmC[ind1] += factor*bbSA;
//				else if (dis.atom2().name().equals("CB"))
//					AtomSumSigmC[ind1] += factor*scSA[dis.atom2().residue().type]/75.0;
//				else
//					AtomSumSigmC[ind1]++;
//				
//				if (dis.atom1().name().equals("CA"))
//					AtomSumSigmC[ind2] += factor*bbSA;
//				else if (dis.atom1().name().equals("CB"))
//					AtomSumSigmC[ind2] += factor*scSA[dis.atom1().residue().type]/75.0;
//				else
//					AtomSumSigmC[ind2]++;
				ind1 = lut[dis.atom1().number()];
				ind2 = lut[dis.atom2().number()];
				AtomSumSigmC[ind1]++;
				AtomSumSigmC[ind2]++;
			}
		}
	}

	//Calculating the energy values. Looping on all the atoms in the protein.
	if (toCalcDerivatives) {
		for (cc=0 ; cc<atomListSize ; cc++) {
			energyVals[cc] = solvCalc.calcRot1(atomTypes[cc], AtomSumSigmC[cc]);
			dEnergyVals[cc] = solvCalc.calcRot1Deriv(atomTypes[cc], AtomSumSigmC[cc]);
			energy += energyVals[cc];
			if (evaluateAtoms)
				atomList.atomAt(cc).addEnergy(energyVals[cc]);
		}
	}
	else {
		for (cc=0 ; cc<atomListSize ; cc++) {
//	Don't need matchingRes now:	if (matchingRes[atomList.atomAt(cc).residueNumber()]){
			energyVals[cc] = solvCalc.calcRot1(atomTypes[cc], AtomSumSigmC[cc]);
			energy += energyVals[cc];
			if (evaluateAtoms)
				atomList.atomAt(cc).addEnergy(energyVals[cc]);
			// Currently, these fields are only assessed without derivatives
			energyHydrophobicWeightedSA += weightsHydrophobicSA[atomTypes[cc]]*energyVals[cc];
			energyPolarWeightedSA += weightsPolarSA[atomTypes[cc]]*energyVals[cc];
			energyRepresentativePolar += weightRepresentativePolars[atomTypes[cc]]*energyVals[cc];
			energyHydrophobicBaysian += weightsHydrophobicMarkerBaysian[atomTypes[cc]]*energyVals[cc];
//			System.out.println(atomList.atomAt(cc).residueNumber()+ " " +
//			 (weightsHydrophobicSA[atomTypes[cc]]*energyVals[cc]+
//					 weightsPolarSA[atomTypes[cc]]*energyVals[cc]));
//	Don't need matchingRes now:		}		
		}
	}
	   	  
	// Second pass over the non-bonded list - for derivatives only
	if (toCalcDerivatives) {
		iter = dislist.iterator();
		while((dis = iter.next()) != null) {
			if (parameters.filterDisForRelevance(dis) && !dis.frozen) {
				sigmaValues = (ROT1SolvationDistanceAttribute) dis.getAttribute(ROT1SolvationDistanceAttribute.SOLVATE_ROT1_ATTRIBUTE);
				ind1 = lut[dis.atom1().number()];
				ind2 = lut[dis.atom2().number()]; 

				// Doing the self derivatives
				forceX[ind1] += dEnergyVals[ind1]*sigmaValues.dsigmCa1dx1;
				forceY[ind1] += dEnergyVals[ind1]*sigmaValues.dsigmCa1dy1;
				forceZ[ind1] += dEnergyVals[ind1]*sigmaValues.dsigmCa1dz1;

				forceX[ind2] += dEnergyVals[ind2]*sigmaValues.dsigmCa2dx2;
				forceY[ind2] += dEnergyVals[ind2]*sigmaValues.dsigmCa2dy2;
				forceZ[ind2] += dEnergyVals[ind2]*sigmaValues.dsigmCa2dz2;

				// Doing the cross derivatives
				forceX[ind2] += dEnergyVals[ind1]*sigmaValues.dsigmCa1dx2;
				forceY[ind2] += dEnergyVals[ind1]*sigmaValues.dsigmCa1dy2;
				forceZ[ind2] += dEnergyVals[ind1]*sigmaValues.dsigmCa1dz2;

				forceX[ind1] += dEnergyVals[ind2]*sigmaValues.dsigmCa2dx1;
				forceY[ind1] += dEnergyVals[ind2]*sigmaValues.dsigmCa2dy1;
				forceZ[ind1] += dEnergyVals[ind2]*sigmaValues.dsigmCa2dz1;
			}
		}  // second pass on the non bonded list                	

		// Finally, the appropriate forces are assigned for every atom. 
		for (cc=0 ; cc<atomListSize ; cc++) {
			atom = atomList.atomAt(cc);
			if (!atom.frozen()) {
				atom.addToFx(-forceX[cc]); // Negating so that it is realy force (and not a mere derivative)
				atom.addToFy(-forceY[cc]); // Negating so that it is realy force
				atom.addToFz(-forceZ[cc]); // Negating so that it is realy force
			}
		}
	}
    return energy;
    }
 
 
    /**
     * Here goes the processing of the various sigmoid functions concerning atom1 and 
     * atom2 in the Distance - dis. The results are updated in the fields of the 
     * SolvateDistanceAttribute of dis - sigmaValues.
     **/ 
    private final void updateSigmVals(Distance dis) {
    	ROT1SolvationDistanceAttribute sigmaValues = 
    		(ROT1SolvationDistanceAttribute) dis.getAttribute(ROT1SolvationDistanceAttribute.SOLVATE_ROT1_ATTRIBUTE); 


    	Sigma.sigma(dis.distance(),end,p1,p2,valAtp1,valAtp2);

    	// Calculating the sigmoid of atom1. 
    	// ---------------------------------------
    	sigmaValues.sigmCa1 = Sigma.s;
    	sigmaValues.dsigmCa1dx1 = Sigma.s_tag*dis.dDistanceDx();
    	sigmaValues.dsigmCa1dy1 = Sigma.s_tag*dis.dDistanceDy();
    	sigmaValues.dsigmCa1dz1 = Sigma.s_tag*dis.dDistanceDz();
    	sigmaValues.dsigmCa1dx2 = -sigmaValues.dsigmCa1dx1;
    	sigmaValues.dsigmCa1dy2 = -sigmaValues.dsigmCa1dy1;
    	sigmaValues.dsigmCa1dz2 = -sigmaValues.dsigmCa1dz1;

    	// Calculating the carbon sigmoid of atom2. 
    	// ---------------------------------------
    	sigmaValues.sigmCa2 = sigmaValues.sigmCa1;
    	sigmaValues.dsigmCa2dx1 = Sigma.s_tag*dis.dDistanceDx();
    	sigmaValues.dsigmCa2dy1 = Sigma.s_tag*dis.dDistanceDy();
    	sigmaValues.dsigmCa2dz1 = Sigma.s_tag*dis.dDistanceDz();
    	sigmaValues.dsigmCa2dx2 = -sigmaValues.dsigmCa2dx1;
    	sigmaValues.dsigmCa2dy2 = -sigmaValues.dsigmCa2dy1;
    	sigmaValues.dsigmCa2dz2 = -sigmaValues.dsigmCa2dz1;

    } // of updateSigmVals    

	public double getEnergyHydrophobicWeightedSA() {
		return energyHydrophobicWeightedSA;
	}

	public double getEnergyPolarWeightedSA() {
		return energyPolarWeightedSA;
	}

	public double getEnergyRepresentativePolar() {
		return energyRepresentativePolar;
	}

	public double getEnergyHydrophobicBaysian() {
		return energyHydrophobicBaysian;
	}
}
	
