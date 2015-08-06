package meshi.energy.solvateNew;
import java.util.Iterator;

import meshi.energy.CooperativeEnergyTerm;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.hydrogenBonds.AbstractHydrogenBond;
import meshi.molecularElements.hydrogenBonds.AbstractHydrogenBondList;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.mathTools.Sigma;
import meshi.util.mathTools.Spline1D;
import meshi.util.mathTools.Spline2D;

/**
 * The implementation of the cooperative solvation term for proteins as described in Kalisman & Keasar (2008). 
 * Since the cooperative solvation is described in the above paper, we bring here only the implementaion
 * details. Especially regarding the calculation of the derivatives, which was too lengthy for the paper. 
 * 
 * The class allows to give a different weight to the solvation energies of certain atom types in the final
 * summation described in Eq. 5. These atom types are side-chain polars, side-chain carbons, and backbone 
 * polars. The class also include a regular hydrogen bond term. The functional form is therefore: 
 *
 * Esolv = weightSCPolarSolvate*Eside_chain_polars + 
 * 		weightSCCarbonSolvate*Eside_chain_carbons + 
 *      weightBBPolarSolvate*Ebackbone_polars + 
 *      weightHB*Ehb
 *      
 * Where Ehb is the negative of the HBC summation over all atoms. The weights are defined in the "Creator" class, 
 * and passed to the constructor as parameters. 
 *
 * General remarks:
 * ----------------
 * 1) This term is derivable twice (because the splines are derivable twice).
 * 2) This term is using hydrogen bond description that is dependent on two angles in the bond. This decription 
 *    follows that of McDonald and Thornton (1994) and Dahiyat et al. (1996). The only place where the hydrogen bond
 *    list is declared explicitly is in line 204. This means that any hydrogen bond implementation that extends the 
 *    "AbstractHydrogenBondList" template can be used, by correcting line 204. 
 * 3) We calculate the regular hydrogen bond energy term (Ehb) together with the solvation terms themselves, since 
 *    the hydrogen bonds calculation is a by-product of the first step in the solvation evaluation, and is thus for free.
 * 4) Disulfide bonds are treated as "hydrogen bonds" between the SG's of two cystines.
 * 5) The SD sulfor of methionine is treated as a hydrophobic carbon.
 * 6) See the remarks in the "Creators" classes for a quick start on how to create a working instance of the term.
 *  
 *
 * The energy evaluation:
 * ----------------------
 * The energy value and derivatives is calculated in 3 steps:
 * 1) A first pass over the non-bonded list. Each Distance instance in the non-bonded-list, is used to update the
 * CNC's and HBC's of its atom pairs (Eqs. 1 and 2, respectively). The partial derivatives of the CNC's and HBC's
 * with respect to the distance atoms are also calculated.  Since some of this values will be needed also in step 3, 
 * we save them in an instance of "SolvateDistanceAttribute" that is attached as an "Attribute" to the Distance instance.
 * 2) A pass on the atom list. Once we have the CNC and HBC of every atom in the protein, we can proceed to calculate the solvation energy 
 * associated with every atom. In this implementation we combined the EI(CNC,HBC) evaluation (Eq. 3) of every atom and 
 * the -log(spline(EI)) evaluation (Eq. 4) into a single step by using a 2D spline, i.e. spline2D(CNC,HBC). The 2D spline 
 * is, of course, atom type specific. The derivatives of each atom solvate energy value with respect to the HBC and CNC  
 * are also calculated.
 * 3) A second pass over the non-bonded list. Equiped with the The derivatives of the atom energies (with respect 
 * to the CNC's and HBC's) from step 2, we can now calculate the energy derivative with respect to the atomic 
 * coordinates. In this step we simply make sure that every term that arises from the derivative chain rule is accounted for. 
 *
 **/
public final class SolvateEnergy extends CooperativeEnergyTerm implements AtomTypes, Residues {
    
	 // Relative strength of SALT BRIDGES compared with regular HYDROGEN BONDS for desolvation purposes, i.e. in
	 // regard to the effect on observed CNC medians. Following Table 1 in the paper. 
	 private final double SALT_BRIDGE_STRENGTH_ASP_OD = 1.0; 	 
	 private final double SALT_BRIDGE_STRENGTH_GLU_OE = 1.0; 	 
	 private final double SALT_BRIDGE_STRENGTH_LYS_NZ = 1.0; 	 
	 private final double SALT_BRIDGE_STRENGTH_ARG_NH = 1.0; 	 
	 private final double SALT_BRIDGE_STRENGTH_TRO = 1.0; 	 
	 private final double SALT_BRIDGE_STRENGTH_TRN = 1.0; 	 
	 /** The following parameter allow for a different weighting of SALT BRIDGES compared with regular HYDROGEN BONDS for the 
	  Ehb energy, that is also claculated. 
	 **/
	 private final double SALT_BRIDGE_STRENGTH_GENERAL = 1.0; 	 

    /** 
     * These are fields for temporary array results that are needed in the evaluation stage.
     * They are declared as fields so that time will not be waisted on creating new 
     * instances of the arrays. The lengths of these arrays is the length of the atom list. 
     * The indexing to these arrays is by the index of the atom in the atom list.
     **/
     private double[] CNC;  // Eq. 1
     private double[] HBC;  // Eq. 2
     private double[] HBCforHBenergy; // May be different because we allow different weighting of the salt bridges in the regular HB energy.
     private double[] dSplineDCNC;     
     private double[] dSplineDHBC;     
     private double[] forceX;
     private double[] forceY;       
     private double[] forceZ;
       
    /** 
     * These following fields are for general use in the class
     **/       
     /** Size of the atom list. **/
     private int atomListSize;  
     /**  The instance of the parameter list object. **/
     private SolvateParametersList parameters; 
     /** The look-up table (lut) converts the atom internal number (field of Atom), which is the index of the array, to its 
      * index in the atom list given to the constructor. **/
     private int[] lut; 	  
     /** Setting the general type for each atom in the atom list: (0) Carbon (1) Backbone polar, (2) Sidechain polar (3) Hydrogens 
      * This is done to save time on type checking. The index to the array is the index of an atom in the atom list. **/
     private int[] superType;  	
     /** The 2D spline array (i.e. spline(CNC,HBC)) for the polar side-chain atoms. The indexing in the array is the 
      * index of the atom type (i.e. 0-189, because we have 190 atom types in meshi currently. The number 190 is not hard coded) **/
     private Spline2D[] splinesSCPolar; 
     /** The 2D spline array (i.e. spline(CNC,HBC)) for the polar backbone atoms. The indexing in the array is the 
      * index of the atom type (i.e. 0-189, because we have 190 atom types in meshi currently. The number 190 is not hard coded) **/
     private Spline2D[] splinesBB; 
     /** The 1D spline array (i.e. spline(CNC,HBC)) for the polar carbon atoms. The indexing in the array is the index of 
      * the atom type (i.e. 0-189, because we have 190 atom types in meshi currently. The number 190 is not hard coded). The splines
      * are 1D because HBC is 0, for carbons. **/
     private Spline1D[] splinesSCCarbon; 
     /** The only hydrogen bond list class in the term. **/ 
     private AbstractHydrogenBondList solvateHB;
     /** Indicating whether the class should also calculate derivatives and forces. These are not necessary, for example, 
      * in side-chain modeling **/ 
     private boolean toCalcDerivatives = true;
     /** We use this array to exclude certain atoms from the solvate calculation, using the method 'inactivateFarFromAtom'. This is
      * useful in side-chain modeling **/ 
     private boolean[] inactive;
     

	 // Weights
	 private double weightSCPolarSolvate;
	 private double weightBBPolarSolvate;
	 private double weightSCCarbonSolvate;
	 private double weightHB;
	 
    public SolvateEnergy() {}
    
    public SolvateEnergy(AtomList atomList, 
            DistanceMatrix dm,
		    SolvateParametersList parameters,
		    AbstractHydrogenBondList hbList,
		    double weightSCPolarSolvate,
		    double weightBBPolarSolvate,
		    double weightSCCarbonSolvate,
		    double weightHB) {
    	this(atomList, dm, parameters, hbList, weightSCPolarSolvate, weightBBPolarSolvate, weightSCCarbonSolvate, weightHB, true);
    }

    
    /** 
     * See the comment at the top of the class for descriptions on the weights.
     **/
    public SolvateEnergy(AtomList atomList, 
                    DistanceMatrix dm,
				    SolvateParametersList parameters,
				    AbstractHydrogenBondList hbList,
				    double weightSCPolarSolvate,
				    double weightBBPolarSolvate,
				    double weightSCCarbonSolvate,
				    double weightHB,
				    boolean toCalcDerivatives) {
	super(toArray(hbList), atomList, dm, parameters, weightSCPolarSolvate);
    this.weightSCPolarSolvate = weightSCPolarSolvate;
    this.weightBBPolarSolvate = weightBBPolarSolvate;
    this.weightSCCarbonSolvate = weightSCCarbonSolvate;
    this.weightHB = weightHB;
    this.toCalcDerivatives = toCalcDerivatives;
    
    int c;
	comment = "Solvation";
	atomListSize = atomList.size();
	this.parameters = parameters;
	if (parameters == null)
	    throw new RuntimeException("The parameters object for this Solvatation term is NULL");
	if (parameters.maxEnd > dm.rMax())
	    throw new RuntimeException("This solvatation term can only work if the rMax in " +
				       "the distance matrix is larger than:" + parameters.maxEnd);
				       
	// Creating the lookup table for the atom numbers.
	// The table converts the atom internal number (field of Atom) to its index  
	// in the atom list given to the constructor.      
	int maxAtomNum = -1;
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


    // Creating the auxilary arrays
    CNC = new double[atomListSize];
    HBC = new double[atomListSize];
    HBCforHBenergy = new double[atomListSize];
    dSplineDCNC = new double[atomListSize];
    dSplineDHBC = new double[atomListSize];
    forceX = new double[atomListSize];
    forceY = new double[atomListSize];
    forceZ = new double[atomListSize];
    superType = new int[atomListSize];
    inactive = new boolean[atomListSize];    

    // All the atoms are active
    for (c=0 ; c<atomListSize ; c++)  
		inactive[c] = false;
		
	// Setting up the splines array, that can associate each residue in the protein with the spline that 
	// suits its type.
	splinesSCPolar = new Spline2D[atomListSize];
	splinesSCCarbon = new Spline1D[atomListSize];
	splinesBB = new Spline2D[atomListSize];
	for (c=0 ; c<atomListSize ; c++)  {
		splinesSCPolar[c] = parameters.scPolarSplines[atomList.atomAt(c).type];
		splinesSCCarbon[c] = parameters.scCarbonSplines[atomList.atomAt(c).type];
		splinesBB[c] = parameters.bbSplines[atomList.atomAt(c).type];
	}
		
	// Setting up the HB list
	solvateHB = hbList;
	
	// Determining the superType for each atom: (0) Hydrophobic, (1) Backbone polar, (2) Sidechain polar (3) Hydrogens
	for (c=0 ; c<atomListSize ; c++)  {
		if (atomList.atomAt(c).isCarbon || (atomList.atomAt(c).type == MSD))
			superType[c] = 0;
		else if ((atomList.atomAt(c).isOxygen || atomList.atomAt(c).isNitrogen) && (atomList.atomAt(c).name().length()==1))
			superType[c] = 1;
		else if ((atomList.atomAt(c).isOxygen || atomList.atomAt(c).isNitrogen || (atomList.atomAt(c).type == CSG)) 
				&& (atomList.atomAt(c).name().length()>1))
			superType[c] = 2;
		else if (atomList.atomAt(c).isHydrogen)
			superType[c] = 3;
		else
			throw new RuntimeException("The following if's should have covered all cases. But due to my carelessness this atom was left out:\n" +
			atomList.atomAt(c) + "\n");
	}
    }  // Of constructor
    
    public void evaluateAtoms() {
		evaluate(true,weightSCPolarSolvate,weightSCCarbonSolvate,weightBBPolarSolvate,weightHB);
    }
    
    /**
     * Calculates Esolv with the weights given in the constructor.
     **/
    public double evaluate() {
    	return evaluate(false,weightSCPolarSolvate,weightSCCarbonSolvate,weightBBPolarSolvate,weightHB);
    }

    /**
     * Calculates Esolv with the weights you give as parameters!
     **/
    public final double evaluate(boolean updateAtoms, double W_SCPolarSolvate, double W_SCCarbonSolvate, double W_BBPolarSolvate,
     							double W_HB) {
    if (! on) return 0.0;
	double energy = 0;
	double energySCPolar=0; 
	double energySCCarbon=0; 
	double energyBBON=0; 
	int cc;        
	DistanceList dislist = dm.nonBondedList();
	Iterator iter;
	Distance dis;
	SolvateDistanceAttribute sigmaValues=null; 
	Atom atom;
	int ind1,ind2;

	//Reseting the auxilary arrays
	for (cc=0 ; cc<atomListSize ; cc++) {
	    CNC[cc] = 0;
	    HBC[cc] = 0; 
	    HBCforHBenergy[cc] = 0; 
	} 
	if (toCalcDerivatives)
		for (cc=0 ; cc<atomListSize ; cc++) {
		    dSplineDCNC[cc] = 0;
		    dSplineDHBC[cc] = 0;
		    forceX[cc] = forceY[cc] = forceZ[cc] = 0.0;
		}

	
	// ***********************************
	// First pass over the non-bonded list
	// ***********************************
	iter = dislist.iterator();
	while((dis = (Distance) iter.next()) != null) {
	   ind1 = lut[dis.atom1().number()];
	   ind2 = lut[dis.atom2().number()];

	   // In case there is not a "SolvateDistanceAttribute" in the distance, we create one.
       if (dis.getAttribute(SolvateDistanceAttribute.SOLVATE_ALL_ATOM_ATTRIBUTE) == null) {
       	  sigmaValues = new SolvateDistanceAttribute();
       	  dis.addAttribute(sigmaValues);
       	  characterizedDistance(dis);
       }
       else
       	  sigmaValues = (SolvateDistanceAttribute) dis.getAttribute(SolvateDistanceAttribute.SOLVATE_ALL_ATOM_ATTRIBUTE);

       if (!inactive[ind1] && !inactive[ind2]) { // Only distances between two active atoms are accounted for
    	   updateSigmVals(dis);
    	   CNC[ind1] += sigmaValues.sigmCa1;
    	   CNC[ind2] += sigmaValues.sigmCa2;
    	   if (sigmaValues.isDisBetween2Polars) {
    		   HBC[ind1] += sigmaValues.saltBridgeFactorA1 * sigmaValues.sigmHB;
    		   HBCforHBenergy[ind1] += sigmaValues.saltBridgeFactorForHBenergyA1 * sigmaValues.sigmHB;
    		   HBC[ind2] += sigmaValues.saltBridgeFactorA2 * sigmaValues.sigmHB;
    		   HBCforHBenergy[ind2] += sigmaValues.saltBridgeFactorForHBenergyA2 * sigmaValues.sigmHB;
    	   }
       }
	}

	// ***************************************************************************
	// A pass over the atom list. Calculating the solvation energies of each atom.
	// ***************************************************************************
	for (cc=0 ; cc<atomListSize ; cc++) 
		if (!inactive[cc]) {
			energySCPolar=energySCCarbon=energyBBON=0.0; 
			if (superType[cc]==0) {		// non-polar atom 
				try {
					splinesSCCarbon[cc].calc(CNC[cc]);
				}
				catch (Exception e) {
					throw new RuntimeException("SCC: " + CNC[cc] + " " + HBC[cc] + "\n" +
							atomList.atomAt(cc) + "\n");
				}
				energySCCarbon = W_SCCarbonSolvate*splinesSCCarbon[cc].s;
				dSplineDCNC[cc] = W_SCCarbonSolvate*splinesSCCarbon[cc].s_tag;
				dSplineDHBC[cc] = 0.0;
			}
			else if (superType[cc]==1) { 		// O and N in the backbone
				try {
					splinesBB[cc].calc(CNC[cc],HBC[cc]);
				}
				catch (Exception e) {
					throw new RuntimeException("BBON: " + CNC[cc] + " " + HBC[cc] + "\n" +
							atomList.atomAt(cc)+ "\n");
				}
				energyBBON = W_BBPolarSolvate*(splinesBB[cc].s);
				dSplineDCNC[cc] = W_BBPolarSolvate*splinesBB[cc].s_tag_x;
				dSplineDHBC[cc] = W_BBPolarSolvate*splinesBB[cc].s_tag_y;
			}
			else if (superType[cc]==2) {	// polar groups in side-chains
				try {
					splinesSCPolar[cc].calc(CNC[cc],HBC[cc]);
				}
				catch (Exception e) {
					throw new RuntimeException("SCP: " + CNC[cc] + " " + HBC[cc] + "\n" +
							atomList.atomAt(cc)+ "\n");
				}
				energySCPolar = W_SCPolarSolvate*splinesSCPolar[cc].s;
				dSplineDCNC[cc] = W_SCPolarSolvate*splinesSCPolar[cc].s_tag_x;
				dSplineDHBC[cc] = W_SCPolarSolvate*splinesSCPolar[cc].s_tag_y;
			}

			energy += (energySCPolar + energySCCarbon + energyBBON - 0.5*W_HB*HBCforHBenergy[cc]); // We want to count each HB once.
			if (updateAtoms) 
				atomList.atomAt(cc).addEnergy(energySCPolar + energySCCarbon + energyBBON - 0.5*W_HB*HBCforHBenergy[cc]);				
		}
	
	if (!toCalcDerivatives)  // From here on its derivatives and forces, so if they are not necessary we can exit now.
		return energy;
   	  
	// ************************************
	// Second pass over the non-bonded list
	// ************************************
	iter = dislist.iterator();
	while((dis = (Distance) iter.next()) != null) {
		ind1 = lut[dis.atom1().number()];
		ind2 = lut[dis.atom2().number()];
		if (!inactive[ind1] && !inactive[ind2]) { // Only distances between two active atoms are accounted for
			sigmaValues = (SolvateDistanceAttribute) dis.getAttribute(SolvateDistanceAttribute.SOLVATE_ALL_ATOM_ATTRIBUTE);
			if (!sigmaValues.isDisInvolvesHydrogens) {
				if (sigmaValues.isDisBetween2Polars && (sigmaValues.sigmHB>0.0)) { // The HB related derivatives 
					solvateHB.findBondByPolars(dis.atom1(),dis.atom2()).applyForcesToAtoms(
							dSplineDHBC[ind1]* sigmaValues.saltBridgeFactorA1 + 
							dSplineDHBC[ind2]* sigmaValues.saltBridgeFactorA2 +
							- 0.5 * W_HB * (sigmaValues.saltBridgeFactorForHBenergyA1 + sigmaValues.saltBridgeFactorForHBenergyA2));
				}
				else {  // Carbon related derivatives 
					// Doing the self derivatives
					forceX[ind1] += dSplineDCNC[ind1] * sigmaValues.dsigmCa1dx1;
					forceY[ind1] += dSplineDCNC[ind1] * sigmaValues.dsigmCa1dy1;
					forceZ[ind1] += dSplineDCNC[ind1] * sigmaValues.dsigmCa1dz1;

					forceX[ind2] += dSplineDCNC[ind2] * sigmaValues.dsigmCa2dx2;
					forceY[ind2] += dSplineDCNC[ind2] * sigmaValues.dsigmCa2dy2;
					forceZ[ind2] += dSplineDCNC[ind2] * sigmaValues.dsigmCa2dz2;

					// Doing the cross derivatives
					forceX[ind2] += dSplineDCNC[ind1] * sigmaValues.dsigmCa1dx2;
					forceY[ind2] += dSplineDCNC[ind1] * sigmaValues.dsigmCa1dy2;
					forceZ[ind2] += dSplineDCNC[ind1] * sigmaValues.dsigmCa1dz2;

					forceX[ind1] += dSplineDCNC[ind2] * sigmaValues.dsigmCa2dx1;
					forceY[ind1] += dSplineDCNC[ind2] * sigmaValues.dsigmCa2dy1;
					forceZ[ind1] += dSplineDCNC[ind2] * sigmaValues.dsigmCa2dz1;   
				}
			}  // No hydrogens are involved	
		}  // Only forces between two activated atoms  
	}   // Of second iteration                 	

    
	// Finally, the appropriate forces are assigned for every atom. 
	for (cc=0 ; cc<atomListSize ; cc++) {
		atom = atomList.atomAt(cc);
		if (!atom.frozen() && !inactive[cc]) {
		    atom.addToFx(-forceX[cc]); // Negating so that it is realy force (and not a mere derivative)
		    atom.addToFy(-forceY[cc]); // Negating so that it is realy force
		    atom.addToFz(-forceZ[cc]); // Negating so that it is realy force
		}
	}
       
    return energy;
    }
 

    /**
     * This method makes changes to the instance of the "SolvateDistanceAttribute" in the distance.
     */
    private final void characterizedDistance(Distance dis) {
    	SolvateDistanceAttribute sigmaValues = (SolvateDistanceAttribute) dis.getAttribute(SolvateDistanceAttribute.SOLVATE_ALL_ATOM_ATTRIBUTE);
    		
    	// Does this distance involve a hydrogen ??
    	// ----------------------------------------
    	if (dis.atom1().isHydrogen || dis.atom2().isHydrogen) 
    		sigmaValues.isDisInvolvesHydrogens = true;
    	else
    		sigmaValues.isDisInvolvesHydrogens = false;
    		
    	// Does this distance involves two polar atoms ??
    	// ----------------------------------------------
    	if ((dis.atom1().isOxygen || dis.atom1().isNitrogen || (dis.atom1().type == CSG)) && 
    		(dis.atom2().isOxygen || dis.atom2().isNitrogen || (dis.atom2().type == CSG))) 
    		sigmaValues.isDisBetween2Polars = true;
    	else 
    		sigmaValues.isDisBetween2Polars = false;

    	if  ((((dis.atom1().type == KNZ) || (dis.atom1().type == RNH) || (dis.atom1().type == TRN)) && 
    	     ((dis.atom2().type == DOD) || (dis.atom2().type == EOE) || (dis.atom2().type == TRO)))      || 
    	     (((dis.atom2().type == KNZ) || (dis.atom2().type == RNH) || (dis.atom2().type == TRN)) && 
    	     ((dis.atom1().type == DOD) || (dis.atom1().type == EOE) || (dis.atom1().type == TRO)))) {  // This is a salt bridge
    	     sigmaValues.saltBridgeFactorForHBenergyA1 = SALT_BRIDGE_STRENGTH_GENERAL;
    	     sigmaValues.saltBridgeFactorForHBenergyA2 = SALT_BRIDGE_STRENGTH_GENERAL;
    	     if (dis.atom1().type == DOD)
    	     	sigmaValues.saltBridgeFactorA1 = SALT_BRIDGE_STRENGTH_ASP_OD;
    	     if (dis.atom1().type == EOE)
    	     	sigmaValues.saltBridgeFactorA1 = SALT_BRIDGE_STRENGTH_GLU_OE;
    	     if (dis.atom1().type == TRO)
    	     	sigmaValues.saltBridgeFactorA1 = SALT_BRIDGE_STRENGTH_TRO;
    	     if (dis.atom1().type == KNZ)
    	     	sigmaValues.saltBridgeFactorA1 = SALT_BRIDGE_STRENGTH_LYS_NZ;
    	     if (dis.atom1().type == RNH)
    	     	sigmaValues.saltBridgeFactorA1 = SALT_BRIDGE_STRENGTH_ARG_NH;    		     
    	     if (dis.atom1().type == TRN)
    	     	sigmaValues.saltBridgeFactorA1= SALT_BRIDGE_STRENGTH_TRN;    		     
    	     if (dis.atom2().type == DOD)
    	     	sigmaValues.saltBridgeFactorA2 = SALT_BRIDGE_STRENGTH_ASP_OD;
    	     if (dis.atom2().type == EOE)
    	     	sigmaValues.saltBridgeFactorA2 = SALT_BRIDGE_STRENGTH_GLU_OE;
    	     if (dis.atom2().type == TRO)
    	     	sigmaValues.saltBridgeFactorA2 = SALT_BRIDGE_STRENGTH_TRO;
    	     if (dis.atom2().type == KNZ)
    	     	sigmaValues.saltBridgeFactorA2 = SALT_BRIDGE_STRENGTH_LYS_NZ;
    	     if (dis.atom2().type == RNH)
    	     	sigmaValues.saltBridgeFactorA2 = SALT_BRIDGE_STRENGTH_ARG_NH;    		     
    	     if (dis.atom2().type == TRN)
    	     	sigmaValues.saltBridgeFactorA2 = SALT_BRIDGE_STRENGTH_TRN;    		     
    	}
    }
 
    /**
     * Here goes the processing of the various sigmoid functions concerning atom1 and 
     * atom2 in the Distance - dis. The results are updated in the fields of the 
     * SolvateDistanceAttribute of dis - sigmaValues.
     **/ 
    private final void updateSigmVals(Distance dis) {
    	SolvateDistanceAttribute sigmaValues = 
    		(SolvateDistanceAttribute) dis.getAttribute(SolvateDistanceAttribute.SOLVATE_ALL_ATOM_ATTRIBUTE); 
       	
    	sigmaValues.resetSigmVals();

    	// Hydrogens are treated as part of the H-bonds
    	// -----------------------------------
    	if (sigmaValues.isDisInvolvesHydrogens) {     	
    		return;
    	}

		// Is this a hydrogen bond or salt bridge?  Possible HB Ahoy!!
    	// -----------------------------------------------------------
		if (sigmaValues.isDisBetween2Polars) {
			AbstractHydrogenBond hb = solvateHB.findBondByPolars(dis.atom1() , dis.atom2());
			if (hb != null) 
				sigmaValues.sigmHB = hb.hbVal();
			else
				sigmaValues.sigmHB = 0.0;
			return;
		}	

    	// Handling CNC
    	// ----------------
		// Atom2 should contribute to atom1's CNC
		Sigma.sigma(dis.distance(),5.5,5.2,5.4,0.95,0.02);
		sigmaValues.sigmCa1 = Sigma.s;
		sigmaValues.dsigmCa1dx1 = Sigma.s_tag*dis.dDistanceDx();
		sigmaValues.dsigmCa1dy1 = Sigma.s_tag*dis.dDistanceDy();
		sigmaValues.dsigmCa1dz1 = Sigma.s_tag*dis.dDistanceDz();
		sigmaValues.dsigmCa1dx2 = -sigmaValues.dsigmCa1dx1;
		sigmaValues.dsigmCa1dy2 = -sigmaValues.dsigmCa1dy1;
		sigmaValues.dsigmCa1dz2 = -sigmaValues.dsigmCa1dz1;
		// Atom1 should contribute to atom2's CNC
		Sigma.sigma(dis.distance(),5.5,5.2,5.4,0.95,0.02);
		sigmaValues.sigmCa2 = Sigma.s;
		sigmaValues.dsigmCa2dx1 = Sigma.s_tag*dis.dDistanceDx();
		sigmaValues.dsigmCa2dy1 = Sigma.s_tag*dis.dDistanceDy();
		sigmaValues.dsigmCa2dz1 = Sigma.s_tag*dis.dDistanceDz();
		sigmaValues.dsigmCa2dx2 = -sigmaValues.dsigmCa2dx1;
		sigmaValues.dsigmCa2dy2 = -sigmaValues.dsigmCa2dy1;
		sigmaValues.dsigmCa2dz2 = -sigmaValues.dsigmCa2dz1;
   	} // of updateSigmVals
    

    public void inactivateFarFromAtom(Atom atom , double R) {
        double x = atom.x();
        double y = atom.y();
        double z = atom.z();
        Atom atom1;
        for (int c=0; c<atomListSize ; c++) {
                atom1 = atomList.atomAt(c);
                if (((atom1.x()-x)>R) || ((atom1.x()-x)<-R) ||
                    ((atom1.y()-y)>R) || ((atom1.y()-y)<-R) ||
                    ((atom1.z()-z)>R) || ((atom1.z()-z)<-R))
                    inactive[c] = true;
                else
                    inactive[c] = false;
        }
    }

}
	

















/*
This is what I put in the spline parameters:

cp scModeling/newForClust/SolvateCarbonSplinesCorrection1.txt ../meshi/parameters/meshiPotential/SolvateCarbonSideChainSplines.dat

cp scModeling/newForClust/SolvatePolarBackboneSplinesOldFashion_noCorr.txt ../meshi/parameters/meshiPotential/SolvatePolarBackboneSplines.dat

cp scModeling/newForClust/SolvatePolarBackboneSplines_NOTBAYES_NoHB.txt ../meshi/parameters/meshiPotential/SolvatePolarBackboneSplines_alt2.dat

cp scModeling/newForClust/SolvatePolarSideChainSplinesOldFashion_noCorr.txt ../meshi/parameters/meshiPotential/SolvatePolarSideChainSplines.dat

cp scModeling/newForClust/SolvatePolarSideChainSplines_NOTBAYES_NoHB.txt ../meshi/parameters/meshiPotential/SolvatePolarSideChainSplines_alt2.dat

*/