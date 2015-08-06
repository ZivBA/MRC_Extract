package meshi.energy.solvateNoHB;
import java.util.Iterator;

import meshi.energy.CooperativeEnergyTerm;
import meshi.energy.CAsolvate.CAsolvateDistanceAttribute;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.parameters.AtomTypes;
import meshi.util.mathTools.Sigma;
import meshi.util.mathTools.Spline1D;


public final class SolvateNoHBEnergy extends CooperativeEnergyTerm implements AtomTypes{

     
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
     private SolvateNoHBParametersList parameters; // The instance of the parameter list object.
     private Spline1D[] splines; // The splines array, that can associate each atom in the protein with the spline that suits its type.
     
    public SolvateNoHBEnergy() {}
    
    public SolvateNoHBEnergy(AtomList atomList, 
                    DistanceMatrix dm,
				    SolvateNoHBParametersList parameters,
				    double weight) {
	super(toArray(dm),atomList, dm, parameters, weight);
	int c;
	int maxAtomNum=-1;
	comment = "Undefined Solvation";
	atomListSize = atomList.size();
	this.parameters = parameters;
	if (parameters == null)
	    throw new RuntimeException("The parameters object for this Solvatation term is NULL");
	if (parameters.maxEnd > dm.rMax())
	    throw new RuntimeException("This solvatation term can only work if the rMax in " +
				       "the distance matrix is larger than:" + parameters.maxEnd);

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
	}
	lut = new int[maxAtomNum+1];
	for (c=0; c<maxAtomNum ; c++) {
	    lut[c] = -1;
	}
	for (c=0; c<atomListSize ; c++) {
	    lut[atomList.atomAt(c).number()] = c;
	}
	
	// Setting up the splines array, that can associate each atom in the protein with the spline that suits its type.
	splines = new Spline1D[atomListSize];
	for (c=0; c<atomListSize ; c++) 
		splines[c] = parameters.atomTypeSplines[atomList.atomAt(c).type];
										
    } // of the constructor
    
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


	//Reseting the auxilary arrays and variables
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
	   envior = AtomSumSigmC[cc]; // The functional form of the environment index
	   splines[cc].calc(envior);
	   atomEnergy = weight*splines[cc].s; // The energy associated with the atom. 
	   dAtomEnergydEnvior[cc] = weight*splines[cc].s_tag; // The derivative of the cooperative part with respect to the environment index.
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
       
    return energy;
    }
 
 
    /**
     * Here goes the processing of the various sigmoid functions concerning atom1 and 
     * atom2 in the Distance - dis. The results are updated in the fields of the 
     * SolvateDistanceAttribute of dis - sigmaValues.
     **/ 
    private final void updateSigmVals(Distance dis) {
    	CAsolvateDistanceAttribute sigmaValues = 
    		(CAsolvateDistanceAttribute) dis.getAttribute(CAsolvateDistanceAttribute.SOLVATE_CA_ATTRIBUTE); 
    	int TsaiAtomicType1 = parameters.atomicTypeConverter[dis.atom1().type]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
    	int TsaiAtomicType2 = parameters.atomicTypeConverter[dis.atom2().type]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
       	
    	sigmaValues.resetAllSigmVals();

    	// Hydrogens are not treated currently
    	// -------------------------------
    	if (dis.atom1().isHydrogen || dis.atom2().isHydrogen) {     	
    		return;
    	}
    	    	

    	// Calculating the carbon sigmoid of atom1. 
    	// ---------------------------------------
    	if (dis.atom2().isCarbon) {
   			Sigma.sigma(dis.distance(),
   						parameters.Cend[TsaiAtomicType1][TsaiAtomicType2],
   						parameters.Cp1[TsaiAtomicType1][TsaiAtomicType2],
   						parameters.Cp2[TsaiAtomicType1][TsaiAtomicType2],
   						parameters.CvalAtp1[TsaiAtomicType1][TsaiAtomicType2],
   						parameters.CvalAtp2[TsaiAtomicType1][TsaiAtomicType2]);
    		sigmaValues.sigmCa1 = Sigma.s;
    		sigmaValues.dsigmCa1dx1 = Sigma.s_tag*dis.dDistanceDx();
  		  	sigmaValues.dsigmCa1dy1 = Sigma.s_tag*dis.dDistanceDy();
  		  	sigmaValues.dsigmCa1dz1 = Sigma.s_tag*dis.dDistanceDz();
 		   	sigmaValues.dsigmCa1dx2 = -sigmaValues.dsigmCa1dx1;
 		   	sigmaValues.dsigmCa1dy2 = -sigmaValues.dsigmCa1dy1;
   		 	sigmaValues.dsigmCa1dz2 = -sigmaValues.dsigmCa1dz1;
   		}

    	// Calculating the carbon sigmoid of atom2. 
    	// ---------------------------------------
    	if (dis.atom1().isCarbon) {	
   			Sigma.sigma(dis.distance(),
   						parameters.Cend[TsaiAtomicType2][TsaiAtomicType1],
   						parameters.Cp1[TsaiAtomicType2][TsaiAtomicType1],
   						parameters.Cp2[TsaiAtomicType2][TsaiAtomicType1],
   						parameters.CvalAtp1[TsaiAtomicType2][TsaiAtomicType1],
   						parameters.CvalAtp2[TsaiAtomicType2][TsaiAtomicType1]);
    		sigmaValues.sigmCa2 = Sigma.s;
    		sigmaValues.dsigmCa2dx1 = Sigma.s_tag*dis.dDistanceDx();
  		  	sigmaValues.dsigmCa2dy1 = Sigma.s_tag*dis.dDistanceDy();
  		  	sigmaValues.dsigmCa2dz1 = Sigma.s_tag*dis.dDistanceDz();
 		   	sigmaValues.dsigmCa2dx2 = -sigmaValues.dsigmCa2dx1;
 		   	sigmaValues.dsigmCa2dy2 = -sigmaValues.dsigmCa2dy1;
   		 	sigmaValues.dsigmCa2dz2 = -sigmaValues.dsigmCa2dz1;
   		}
    	
   	} // of updateSigmVals
    
}
	
