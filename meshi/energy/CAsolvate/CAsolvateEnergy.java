package meshi.energy.CAsolvate;
import java.util.Iterator;

import meshi.energy.CooperativeEnergyTerm;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.mathTools.Sigma;

/**
 * The CAsolvateEnergy gives high (but derivable) energy values if in a radius 'Rad' around a certain CA atom
 * there are more than 'MaxCANeighbors' CA atoms. If the number of neighbors does not exceed 'MaxCANeighbors' 
 * the energy is 0.
 *
 * General remarks:
 * ----------------
 * 1) This term is derivable once.
 * 2) Here we bring only the implementaion details. Especially regarding the calculation of the derivatives.
 * 3) See the remark on the class constructor, and the remarks in CAsolvateCreator for a quick start on how 
 * to create a working instance of the term.
 *  
 *
 * The energy evaluation:
 * ----------------------
 * The energy value and derivatives is calculated in 3 steps:
 * 1) A first pass over the non-bonded list. For each Distance in the list, the carbon sigmoid 
 * values are calculated along with their derivatives. (Since some of this values will be needed also in step 3, 
 * we save them in an instance of CAsolvateDistanceAttribute that is attached as an attribute to the Distance 
 * instance.) The sigmoid values are used to update the carbon neighbor count of the two atoms in the Distance.
 * 2) Once we have the carbon neighbor count we can proceed to calculate the energy assiciate with each atom. If 
 * the number of neighbours is less than 'MaxCANeighbors' the energy is 0. If the number of neighbors is higher, the
 * energy goes up quadratically in (number of neighbors - MaxCANeighbors). The derivatives of the energy with 
 * respect to the number of carcon neighbors are also calculated.
 * 3) A second pass over the non-bonded list. Equiped with the The derivatives of the atom energies (with respect 
 * to the number of carcon neighbors) from step 2, we can now calculate the energy derivative with respect to the atomic 
 * coordinates. Note, that there are 2 types of derivatives. "Self derivatives" arises from the affect of an atom 
 * coordinates on its own sigmoid values. "Cross derivatives" arises from the affect of an atom coordinates on 
 * the sigmoid values of its partner in the Distance. 
 *
 **/
public final class CAsolvateEnergy extends CooperativeEnergyTerm {

     
    /** 
     * These are fields for temporary array results that are needed in the evaluation stage.
     * They are declared as fields so that time will not be waisted on creating new 
     * instances of the arrays.
     **/
     private double[] AtomSumSigmC;
     private double[] dAtomEnergydEnvior;
     private double[] forceX;
     private double[] forceY;       
     private double[] forceZ;
       
    /** 
     * These fields are for general use in the class
     **/       
     private int[] lut; 	// The look-up table (lut) converts the atom internal number (field of Atom) to its index in the atom list given to the constructor.      
     private int atomListSize;
     private double Rad; 
     private double MaxCANeighbors; 
     private double cEnd,p1,p2,valAtp1=0.98,valAtp2=0.02;
     
    public CAsolvateEnergy() {}
    
/**
 * The constructor parameters:
 * atomList,dm,weight - standard CooperativeEnergyTerm inputs. 
 * simpleHBweight - The weight of the HB part.
 * sigmoidBeginsWithH,sigmoidEndsWithH - The transition angles (in degrees) for the HB sigmoid when one of the 
 * base atoms (or both) is a hydrogen. Above sigmoidEndsWithH the sigmoid is given a value of 1.0 . Bellow          
 * sigmoidBeginsWithH the sigmoid is given a value of 0.0 . In between it raises smoothly by cubic spline.
 * sigmoidBeginsNoH,sigmoidEndsNoH - The same as above, only for HB sigmoids where none of the base atoms
 * is a hydrogen. 
 **/
    public CAsolvateEnergy(AtomList atomList, 
                    DistanceMatrix dm,
				    CAsolvateParametersList parameters,
				    double Rad,
			    	double MaxCANeighbors, 
				    double weight) {
	super(toArray(dm),atomList, dm, parameters, weight);
	int c;
	int maxAtomNum=-1;
	comment = "Undefined Solvation";
	atomListSize = atomList.size();
	this.Rad = Rad;
	this.MaxCANeighbors = MaxCANeighbors;
	if (Rad*1.05 > dm.rMax())
	    throw new RuntimeException("This solvatation term can only work if the rMax in " +
				       "the distance matrix is larger than (Rad*1.05):" + (Rad*1.05));
	cEnd = dm.rMax();
	p1 = Rad*0.95;
	p2 = Rad*1.05;

    // Creating the auxilary arrays
    AtomSumSigmC = new double[atomListSize];
    dAtomEnergydEnvior = new double[atomListSize];
    forceX = new double[atomListSize];
    forceY = new double[atomListSize];
    forceZ = new double[atomListSize];
    
    
	// Creating the lookup table for the atom numbers.
	// The table converts the atom internal number (field of Atom) to its index  
	// in the atom list given to the constructor.      
	for (c=0; c<atomListSize ; c++) {
	    if (atomList.atomAt(c).number() > maxAtomNum)
		maxAtomNum = atomList.atomAt(c).number();
		if (!atomList.atomAt(c).name().equals("CA"))
		   throw new RuntimeException("\n\nThe atom list to the CA solave must be composed of CA atoms only.\n" +
		   "However, The following atom found:\n" +  atomList.atomAt(c) + "\n\n");
	}
	lut = new int[maxAtomNum+1];
	for (c=0; c<maxAtomNum ; c++) {
	    lut[c] = -1;
	}
	for (c=0; c<atomListSize ; c++) {
	    lut[atomList.atomAt(c).number()] = c;
	}	
    }
    
    public void setComment(String str) {
    	comment = str;
    }
     
    public void evaluateAtoms() {
		evaluate(true);
    }
    
    public double evaluate() {
    	return evaluate(false);
    }

    public final double evaluate(boolean updateAtoms) {
	if (! on) return 0.0;
	double energy = 0;
	double envior;
	double atomEnergy=0; 
	int cc;        
	DistanceList dislist = dm.nonBondedList();
	Iterator iter;
	Distance dis;
	CAsolvateDistanceAttribute sigmaValues; 
	Atom atom;
	int ind1,ind2;


	//Reseting the auxilary arrays
	for (cc=0 ; cc<atomListSize ; cc++) {
	    AtomSumSigmC[cc] = 0; 
	    forceX[cc] = forceY[cc] = forceZ[cc] = 0.0;
	} 
	
	// First pass over the non-bonded list
	iter = dislist.iterator();
	while((dis = (Distance) iter.next()) != null) {
       if (dis.getAttribute(CAsolvateDistanceAttribute.SOLVATE_CA_ATTRIBUTE) == null) {
       	  sigmaValues = new CAsolvateDistanceAttribute();
       	  dis.addAttribute(sigmaValues);
       }
       else
          sigmaValues = (CAsolvateDistanceAttribute) dis.getAttribute(CAsolvateDistanceAttribute.SOLVATE_CA_ATTRIBUTE);
	   ind1 = lut[dis.atom1().number()];
	   ind2 = lut[dis.atom2().number()];
	   updateSigmVals(dis);
	   AtomSumSigmC[ind1] += sigmaValues.sigmCa1;
	   AtomSumSigmC[ind2] += sigmaValues.sigmCa2;
	}

	//Calculating the energy values. Looping on all the atoms in the protein.
	for (cc=0 ; cc<atomListSize ; cc++) {
		if (AtomSumSigmC[cc]<MaxCANeighbors) {
			atomEnergy = 0.0;
			dAtomEnergydEnvior[cc] = 0.0;
		}
		else {
			atomEnergy = weight*(AtomSumSigmC[cc]-MaxCANeighbors)*(AtomSumSigmC[cc]-MaxCANeighbors);
			dAtomEnergydEnvior[cc] = weight*2.0*(AtomSumSigmC[cc]-MaxCANeighbors);			
		}		
	   energy += atomEnergy;
	   if (updateAtoms) 
	      	atomList.atomAt(cc).addEnergy(atomEnergy);
	}
	   	  
	// Second pass over the non-bonded list
	iter = dislist.iterator();
	while((dis = (Distance) iter.next()) != null) {
       sigmaValues = (CAsolvateDistanceAttribute) dis.getAttribute(CAsolvateDistanceAttribute.SOLVATE_CA_ATTRIBUTE);
	   ind1 = lut[dis.atom1().number()];
	   ind2 = lut[dis.atom2().number()]; 
	   // Doing the self derivatives
	   forceX[ind1] += dAtomEnergydEnvior[ind1]*sigmaValues.dsigmCa1dx1;
	   forceY[ind1] += dAtomEnergydEnvior[ind1]*sigmaValues.dsigmCa1dy1;
	   forceZ[ind1] += dAtomEnergydEnvior[ind1]*sigmaValues.dsigmCa1dz1;
	   
	   forceX[ind2] += dAtomEnergydEnvior[ind2]*sigmaValues.dsigmCa2dx2;
	   forceY[ind2] += dAtomEnergydEnvior[ind2]*sigmaValues.dsigmCa2dy2;
	   forceZ[ind2] += dAtomEnergydEnvior[ind2]*sigmaValues.dsigmCa2dz2; 
	   
	   // Doing the cross derivatives
	   forceX[ind2] += dAtomEnergydEnvior[ind1]*sigmaValues.dsigmCa1dx2;
	   forceY[ind2] += dAtomEnergydEnvior[ind1]*sigmaValues.dsigmCa1dy2;
	   forceZ[ind2] += dAtomEnergydEnvior[ind1]*sigmaValues.dsigmCa1dz2;
	                  
	   forceX[ind1] += dAtomEnergydEnvior[ind2]*sigmaValues.dsigmCa2dx1; 
	   forceY[ind1] += dAtomEnergydEnvior[ind2]*sigmaValues.dsigmCa2dy1; 
	   forceZ[ind1] += dAtomEnergydEnvior[ind2]*sigmaValues.dsigmCa2dz1;
	}                 	
    
	// Finally, the appropriate forces are assigned for every atom. 
	for (cc=0 ; cc<atomListSize ; cc++) {
		atom = atomList.atomAt(cc);
		if (! atom.frozen()) {
		    atom.addToFx(-forceX[cc]); // Negating so that it is realy force (and not a mere derivative)
		    atom.addToFy(-forceY[cc]); // Negating so that it is realy force
		    atom.addToFz(-forceZ[cc]); // Negating so that it is realy force
		}
	}
       
    return energy;
    }
 
 
    /**
     * Here goes the processing of the various sigmoid functions concerning atom1 and 
     * atom2 in the Distance - dis. The results are updated in the fields of the 
     * CAsolvateDistanceAttribute of dis.
     **/ 
    private final void updateSigmVals(Distance dis) {
    	CAsolvateDistanceAttribute sigmaValues = 
    		(CAsolvateDistanceAttribute) dis.getAttribute(CAsolvateDistanceAttribute.SOLVATE_CA_ATTRIBUTE); 
       	
    	// Calculating the carbon sigmoids 
    	// ---------------------------------------
   			Sigma.sigma(dis.distance(),cEnd,p1,p2,valAtp1,valAtp2);
    		sigmaValues.sigmCa1 = Sigma.s;
    		sigmaValues.dsigmCa1dx1 = Sigma.s_tag*dis.dDistanceDx();
  		  	sigmaValues.dsigmCa1dy1 = Sigma.s_tag*dis.dDistanceDy();
  		  	sigmaValues.dsigmCa1dz1 = Sigma.s_tag*dis.dDistanceDz();
 		   	sigmaValues.dsigmCa1dx2 = -sigmaValues.dsigmCa1dx1;
 		   	sigmaValues.dsigmCa1dy2 = -sigmaValues.dsigmCa1dy1;
   		 	sigmaValues.dsigmCa1dz2 = -sigmaValues.dsigmCa1dz1;
    		sigmaValues.sigmCa2 = sigmaValues.sigmCa1;
    		sigmaValues.dsigmCa2dx1 = sigmaValues.dsigmCa1dx1;
  		  	sigmaValues.dsigmCa2dy1 = sigmaValues.dsigmCa1dy1;
  		  	sigmaValues.dsigmCa2dz1 = sigmaValues.dsigmCa1dz1;
 		   	sigmaValues.dsigmCa2dx2 = sigmaValues.dsigmCa1dx2;
 		   	sigmaValues.dsigmCa2dy2 = sigmaValues.dsigmCa1dy2;
   		 	sigmaValues.dsigmCa2dz2 = sigmaValues.dsigmCa1dz2;

   	} // of updateSigmVals


    
}
	
