package meshi.energy.oldSolvate;
import java.util.Iterator;

import meshi.energy.CooperativeEnergyTerm;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.parameters.AtomTypes;
import meshi.util.mathTools.Sigma;
import meshi.util.mathTools.Spline1D;

/**
 * The implementation of the cooperative solvation term for proteins described in Kalisman & Keasar (2006). 
 * The functional form is:
 *
 * Eterm = Wcooperative*Ecoop - Whb*Ehb
 *
 * General remarks:
 * ----------------
 * 1) This term is derivable once.
 * 2) Since the cooperative solvation is described in the above paper, we bring here only the implementaion
 * details. Especially regarding the calculation of the derivatives, which was too lengthy for the paper. 
 * 3) We calculate the regular hydrogen bond energy term (Ehb) together with the cooperative term (Ecoop) 
 * since the hydrogen bonds calculation is a by-product of the first step in the Ecoop evaluation.
 * 4) The hydrogen bonds are angle dependent. See the remark on SolvateHBAngle for the definintion of these
 * angles.
 * 5) Disulfide bonds are treated as "hydrogen bonds" between the CG's of two cystines.
 * 6) The weight of the cooperative part (Wcooperative) is the standard MESHI weight. The HB weight (Whb) is 
 * set either as a parameter to the constructor, or later in a designated method: 'setSimpleHBweight'
 * 7) See the remark on the class constructor, and the remarks in SolvateCreator for a quick start on how 
 * to create a working instance of the term.
 *  
 *
 * The energy evaluation:
 * ----------------------
 * The energy value and derivatives is calculated in 3 steps:
 * 1) A first pass over the non-bonded list. For each Distance in the list, the carbon and HB sigmoid 
 * values are calculated along with their derivatives. (Since some of this values will be needed also in step 3, 
 * we save them in an instance of SolvateDistanceAttribute that is attached as an attribute to the Distance 
 * instance.) The sigmoid values are used to update the carbon neighbor count and the HB count of the two atoms    
 * in the Distance.
 * 2) Once we have the carbon neighbor count and the HB count of every atom in the protein, we can proceed to 
 * calculate the solvation energy associated with every atom. By looping on all the atoms in the protein, the 
 * environment index of each atom is calculated. The environment index value is passed into a 1D spline that 
 * turns it into an energy value. The spline is, of course, atom type specific. The derivative of the atom 
 * energy value with respect to the environment index are also calculated.
 * 3) A second pass over the non-bonded list. Equiped with the The derivatives of the atom energies (with respect 
 * to the environment indices) from step 2, we can now calculate the energy derivative with respect to the atomic 
 * coordinates. Note, that there are 3 types of derivatives. "Self derivatives" arises from the affect of an atom 
 * coordinates on its own sigmoid values. "Cross derivatives" arises from the affect of an atom coordinates on 
 * the sigmoid values of its partner in the Distance. "HB angle derivatives" arises from the affect the base atom 
 * coordinates on the HB sigmoids of polar atoms participating in hydrogen bonding.  
 *
 **/
public final class SolvateEnergy extends CooperativeEnergyTerm implements AtomTypes{
private boolean updateA;
    /** 
     * These are fields for temporary array results that are needed in the evaluation stage.
     * They are declared as fields so that time will not be waisted on creating new 
     * instances of the arrays.
     **/
     private double[] AtomSumSigmC;
     private double[] AtomSumSigmHB;
     private double[] AtomSumSigmPHB;
     private double[] dAtomEnergydEnvior;
     private double[] forceX;
     private double[] forceY;       
     private double[] forceZ;
       
    /** 
     * These fields are for general use in the class
     **/       
     private int[] lut; 	// The look-up table (lut) converts the atom internal number (field of Atom) to its index in the atom list given to the constructor.      
     private int atomListSize;
     private SolvateParametersList parameters; // The instance of the parameter list object.
     private double simpleHBweight; // The weight given to the HB energy term.
     private Spline1D[] splines; // The splines array, that can associate each atom in the protein with the spline that suits its type.
     private SolvateHBAngle solvateHBAngle; // We create a single instance of SolvateHBAngle, and use it to calculate values, regarding the angular properties of the hydrogen bonds.
     private Atom[] baseAtom; // This array stores pointers to the base atoms of every non-hydrogen polar atom in the protein.   
     
    public SolvateEnergy() {}
    
/**
 * The constructor parameters:
 * atomList,dm,weight - standard CooperativeEnergyTerm inputs. The 'weight' parameter is the parmeter of 
 * the cooperative part (Wcooperative). Note, that changing it will not affect the weight given to the HB 
 * part.
 * simpleHBweight - The weight of the HB part.
 * sigmoidBeginsWithH,sigmoidEndsWithH - The transition angles (in degrees) for the HB sigmoid when one of the 
 * base atoms (or both) is a hydrogen. Above sigmoidEndsWithH the sigmoid is given a value of 1.0 . Bellow          
 * sigmoidBeginsWithH the sigmoid is given a value of 0.0 . In between it raises smoothly by cubic spline.
 * sigmoidBeginsNoH,sigmoidEndsNoH - The same as above, only for HB sigmoids where none of the base atoms
 * is a hydrogen. 
 **/
    public SolvateEnergy(AtomList atomList, 
                    DistanceMatrix dm,
				    SolvateParametersList parameters,
				    double simpleHBweight,
			    	double sigmoidBeginsWithH, 
    				double sigmoidEndsWithH,
   				 	double sigmoidBeginsNoH,
    				double sigmoidEndsNoH, 
				    double weight) {
	super(toArray(),atomList, dm, parameters, weight);
	int c;
	int maxAtomNum=-1;
	comment = "Undefined Solvation";
	atomListSize = atomList.size();
	this.simpleHBweight = simpleHBweight;
	this.parameters = parameters;
	if (parameters == null)
	    throw new RuntimeException("The parameters object for this Solvatation term is NULL");
	if (parameters.maxEnd > dm.rMax())
	    throw new RuntimeException("This solvatation term can only work if the rMax in " +
				       "the distance matrix is larger than:" + parameters.maxEnd);

    // Creating the auxilary arrays
    AtomSumSigmC = new double[atomListSize];
    AtomSumSigmHB = new double[atomListSize];
    AtomSumSigmPHB = new double[atomListSize];
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
		
	// Setting up the HB angle object, and the base atom array. 
	// See the documentaion of "SolvateHBAngle.java" to see what are the parameters to its
	// constructor. 
	solvateHBAngle = new SolvateHBAngle(dm,
						sigmoidBeginsWithH,
						sigmoidEndsWithH,
						sigmoidBeginsNoH,
						sigmoidEndsNoH);
							
	// Setting the base atom for every polar non-hydrogen atom. This base atom is used to determined the HB angles.
	// The places in 'baseAtom' that correspond to carbon or hydrogen atoms remain 'null'.  
	baseAtom = new Atom[atomListSize];
	setBaseAtom();
    }
    
    public void setComment(String str) {
    	comment = str;
    }
    
    public void setSimpleHBweight(double newWeight){
    	simpleHBweight = newWeight;
    }
 
    public void update(int numberOfUpdates) {}

    public void evaluateAtoms() {
		evaluate(true,weight,simpleHBweight);
    }
    
    public double evaluate() {
    	return evaluate(false,weight,simpleHBweight);
    }

    public final double evaluate(boolean updateAtoms, double cooperativeWeight, double simpleHBweight) {
//	updateA = updateAtoms;
        if (! on) return 0.0;
	double energy = 0;
	double envior;
	double atomEnergy=0; 
	int cc;        
	DistanceList dislist = dm.nonBondedList();
	Iterator iter;
	Distance dis;
	SolvateDistanceAttribute sigmaValues; 
	Atom atom;
	int ind1,ind2,ind3,ind4;


	//Reseting the auxilary arrays
	for (cc=0 ; cc<atomListSize ; cc++) {
	    AtomSumSigmC[cc] = 0; 
	    AtomSumSigmHB[cc] = 0; 
	    AtomSumSigmPHB[cc] = 0; 
	    forceX[cc] = forceY[cc] = forceZ[cc] = 0.0;
	} 
	
	// First pass over the non-bonded list
	iter = dislist.iterator();
	while((dis = (Distance) iter.next()) != null) {
       if (dis.getAttribute(SolvateDistanceAttribute.SOLVATE_ALL_ATOM_ATTRIBUTE) == null) {
       	  sigmaValues = new SolvateDistanceAttribute();
       	  dis.addAttribute(sigmaValues);
       }
       else
          sigmaValues = (SolvateDistanceAttribute) dis.getAttribute(SolvateDistanceAttribute.SOLVATE_ALL_ATOM_ATTRIBUTE);
	   ind1 = lut[dis.atom1().number()];
	   ind2 = lut[dis.atom2().number()];
	   updateSigmVals(dis);
	   AtomSumSigmC[ind1] += sigmaValues.sigmCa1;
	   AtomSumSigmC[ind2] += sigmaValues.sigmCa2;
	   AtomSumSigmHB[ind1] += sigmaValues.sigmHBa1;
	   AtomSumSigmHB[ind2] += sigmaValues.sigmHBa2;
	   AtomSumSigmPHB[ind1] += sigmaValues.sigmPHBa1;
	   AtomSumSigmPHB[ind2] += sigmaValues.sigmPHBa2;
/*
if (updateAtoms && ((dis.atom1().residueNumber()==155) && dis.atom1().name().equals("OE1") ||
(dis.atom2().residueNumber()==155) && dis.atom2().name().equals("OE1"))) {
dis.atoms().print();
System.out.println("\n" +sigmaValues.sigmHBa1+ " " + sigmaValues.sigmHBa2 + " " +  sigmaValues.sigmPHBa1 
+ " " +  sigmaValues.sigmPHBa2 + "\n");
}
if (updateAtoms && ((dis.atom1().residueNumber()==155) && dis.atom1().name().equals("OE2") ||
(dis.atom2().residueNumber()==155) && dis.atom2().name().equals("OE2"))) {
dis.atoms().print();
System.out.println("\n" +sigmaValues.sigmHBa1+ " " + sigmaValues.sigmHBa2 + " " +  sigmaValues.sigmPHBa1
+ " " +  sigmaValues.sigmPHBa2 + "\n");
}
*/
	}

	//Calculating the energy values. Looping on all the atoms in the protein.
	for (cc=0 ; cc<atomListSize ; cc++) {
	   envior = AtomSumSigmC[cc]*(1-AtomSumSigmHB[cc]); // The functional form of the environment index
/*
if (updateAtoms && (atomList.atomAt(cc).residueNumber()==155) && atomList.atomAt(cc).name().equals("OE1"))
System.out.println(atomList.atomAt(cc) + "\n" +  envior + " " +  AtomSumSigmC[cc] + " " + AtomSumSigmHB[cc]);
if (updateAtoms && (atomList.atomAt(cc).residueNumber()==155) && atomList.atomAt(cc).name().equals("OE2"))
System.out.println(atomList.atomAt(cc) + "\n" +  envior + " " +  AtomSumSigmC[cc] + " " + AtomSumSigmHB[cc]);
*/
	   if (envior<0.0)
	      envior = 0.0;
	   splines[cc].calc(envior);
	   atomEnergy = cooperativeWeight*splines[cc].s - simpleHBweight*AtomSumSigmPHB[cc]; // The energy associated with the atom. 
	   dAtomEnergydEnvior[cc] = cooperativeWeight*splines[cc].s_tag; // The derivative of the cooperative part with respect to the environment index.
	   energy += atomEnergy;
	   if (updateAtoms) 
	      	atomList.atomAt(cc).addEnergy(atomEnergy);
	}
	   	  
	// Second pass over the non-bonded list
	iter = dislist.iterator();
	while((dis = (Distance) iter.next()) != null) {
       sigmaValues = (SolvateDistanceAttribute) dis.getAttribute(SolvateDistanceAttribute.SOLVATE_ALL_ATOM_ATTRIBUTE);
	   ind1 = lut[dis.atom1().number()];
	   ind2 = lut[dis.atom2().number()]; 
	   // Doing the self derivatives
	   forceX[ind1] += dAtomEnergydEnvior[ind1]*
	                 ((1-AtomSumSigmHB[ind1])*sigmaValues.dsigmCa1dx1 -
	                  AtomSumSigmC[ind1]*sigmaValues.dsigmHBa1dx1) - 
	                  simpleHBweight*sigmaValues.dsigmPHBa1dx1;
	   forceY[ind1] += dAtomEnergydEnvior[ind1]*
	                 ((1-AtomSumSigmHB[ind1])*sigmaValues.dsigmCa1dy1 -
	                  AtomSumSigmC[ind1]*sigmaValues.dsigmHBa1dy1) - 
	                  simpleHBweight*sigmaValues.dsigmPHBa1dy1;
	   forceZ[ind1] += dAtomEnergydEnvior[ind1]*
	                 ((1-AtomSumSigmHB[ind1])*sigmaValues.dsigmCa1dz1 -
	                  AtomSumSigmC[ind1]*sigmaValues.dsigmHBa1dz1) - 
	                  simpleHBweight*sigmaValues.dsigmPHBa1dz1;
	   
	   forceX[ind2] += dAtomEnergydEnvior[ind2]*
	                 ((1-AtomSumSigmHB[ind2])*sigmaValues.dsigmCa2dx2 -
	                  AtomSumSigmC[ind2]*sigmaValues.dsigmHBa2dx2) - 
	                  simpleHBweight*sigmaValues.dsigmPHBa2dx2;
	   forceY[ind2] += dAtomEnergydEnvior[ind2]*
	                 ((1-AtomSumSigmHB[ind2])*sigmaValues.dsigmCa2dy2 -
	                  AtomSumSigmC[ind2]*sigmaValues.dsigmHBa2dy2) - 
	                  simpleHBweight*sigmaValues.dsigmPHBa2dy2;	
	   forceZ[ind2] += dAtomEnergydEnvior[ind2]*
	                 ((1-AtomSumSigmHB[ind2])*sigmaValues.dsigmCa2dz2 -
	                  AtomSumSigmC[ind2]*sigmaValues.dsigmHBa2dz2) - 
	                  simpleHBweight*sigmaValues.dsigmPHBa2dz2;

	   // Doing the cross derivatives
	   forceX[ind2] += dAtomEnergydEnvior[ind1]*
	                 ((1-AtomSumSigmHB[ind1])*sigmaValues.dsigmCa1dx2 -
	                  AtomSumSigmC[ind1]*sigmaValues.dsigmHBa1dx2) - 
	                  simpleHBweight*sigmaValues.dsigmPHBa1dx2;
	   forceY[ind2] += dAtomEnergydEnvior[ind1]*
	                 ((1-AtomSumSigmHB[ind1])*sigmaValues.dsigmCa1dy2 -
	                  AtomSumSigmC[ind1]*sigmaValues.dsigmHBa1dy2) - 
	                  simpleHBweight*sigmaValues.dsigmPHBa1dy2;
	   forceZ[ind2] += dAtomEnergydEnvior[ind1]*
	                 ((1-AtomSumSigmHB[ind1])*sigmaValues.dsigmCa1dz2 -
	                  AtomSumSigmC[ind1]*sigmaValues.dsigmHBa1dz2) - 
	                  simpleHBweight*sigmaValues.dsigmPHBa1dz2;
	                  
	   forceX[ind1] += dAtomEnergydEnvior[ind2]*
	                 ((1-AtomSumSigmHB[ind2])*sigmaValues.dsigmCa2dx1 -
	                  AtomSumSigmC[ind2]*sigmaValues.dsigmHBa2dx1) - 
	                  simpleHBweight*sigmaValues.dsigmPHBa2dx1;
	   forceY[ind1] += dAtomEnergydEnvior[ind2]*
	                 ((1-AtomSumSigmHB[ind2])*sigmaValues.dsigmCa2dy1 -
	                  AtomSumSigmC[ind2]*sigmaValues.dsigmHBa2dy1) - 
	                  simpleHBweight*sigmaValues.dsigmPHBa2dy1;	
	   forceZ[ind1] += dAtomEnergydEnvior[ind2]*
	                 ((1-AtomSumSigmHB[ind2])*sigmaValues.dsigmCa2dz1 -
	                  AtomSumSigmC[ind2]*sigmaValues.dsigmHBa2dz1) - 
	                  simpleHBweight*sigmaValues.dsigmPHBa2dz1;
	                  
	   // Doing the derivatives that arises from the angular part of the HB involving 4 atoms
	   if (sigmaValues.hbAngleDerivativeNonZero) {	
	   		ind3 = lut[baseAtom[lut[dis.atom1().number()]].number()];
	   		ind4 = lut[baseAtom[lut[dis.atom2().number()]].number()];
	   			   		
			forceX[ind3] -= (dAtomEnergydEnvior[ind1]*AtomSumSigmC[ind1]*sigmaValues.dsigmHBa1dx3 +
		   					dAtomEnergydEnvior[ind2]*AtomSumSigmC[ind2]*sigmaValues.dsigmHBa2dx3 +
		   					simpleHBweight*(sigmaValues.dsigmPHBa1dx3+sigmaValues.dsigmPHBa2dx3));
	   		
			forceY[ind3] -= (dAtomEnergydEnvior[ind1]*AtomSumSigmC[ind1]*sigmaValues.dsigmHBa1dy3 +
		   					dAtomEnergydEnvior[ind2]*AtomSumSigmC[ind2]*sigmaValues.dsigmHBa2dy3 +
		   					simpleHBweight*(sigmaValues.dsigmPHBa1dy3+sigmaValues.dsigmPHBa2dy3));
	   		
			forceZ[ind3] -= (dAtomEnergydEnvior[ind1]*AtomSumSigmC[ind1]*sigmaValues.dsigmHBa1dz3 +
		   					dAtomEnergydEnvior[ind2]*AtomSumSigmC[ind2]*sigmaValues.dsigmHBa2dz3 +
		   					simpleHBweight*(sigmaValues.dsigmPHBa1dz3+sigmaValues.dsigmPHBa2dz3));
	   		
			forceX[ind4] -= (dAtomEnergydEnvior[ind1]*AtomSumSigmC[ind1]*sigmaValues.dsigmHBa1dx4 +
		   					dAtomEnergydEnvior[ind2]*AtomSumSigmC[ind2]*sigmaValues.dsigmHBa2dx4 +
		   					simpleHBweight*(sigmaValues.dsigmPHBa1dx4+sigmaValues.dsigmPHBa2dx4));
	   		
			forceY[ind4] -= (dAtomEnergydEnvior[ind1]*AtomSumSigmC[ind1]*sigmaValues.dsigmHBa1dy4 +
		   					dAtomEnergydEnvior[ind2]*AtomSumSigmC[ind2]*sigmaValues.dsigmHBa2dy4 +
		   					simpleHBweight*(sigmaValues.dsigmPHBa1dy4+sigmaValues.dsigmPHBa2dy4));
	   		
			forceZ[ind4] -= (dAtomEnergydEnvior[ind1]*AtomSumSigmC[ind1]*sigmaValues.dsigmHBa1dz4 +
		   					dAtomEnergydEnvior[ind2]*AtomSumSigmC[ind2]*sigmaValues.dsigmHBa2dz4 +
		   					simpleHBweight*(sigmaValues.dsigmPHBa1dz4+sigmaValues.dsigmPHBa2dz4));				
	   }
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
     * SolvateDistanceAttribute of dis - sigmaValues.
     **/ 
    private final void updateSigmVals(Distance dis) {
    	SolvateDistanceAttribute sigmaValues = 
    		(SolvateDistanceAttribute) dis.getAttribute(SolvateDistanceAttribute.SOLVATE_ALL_ATOM_ATTRIBUTE); 
    	int TsaiAtomicType1 = parameters.atomicTypeConverter[dis.atom1().type]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
    	int TsaiAtomicType2 = parameters.atomicTypeConverter[dis.atom2().type]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
    	double tmpSigmaHBatom1,tmpSigma_TagHBatom1,tmpSigmaHBatom2,tmpSigma_TagHBatom2;
    	double tmpSigmaPHBatom1,tmpSigma_TagPHBatom1,tmpSigmaPHBatom2,tmpSigma_TagPHBatom2;
       	
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
    	
    	// Possible HB Ahoy!!
    	// (atom1 and atom2 are capable of forming a hydrogen bond) 
    	// ------------------
    	if ((dis.atom1().isOxygen || dis.atom1().isNitrogen || dis.atom1().isSulfur) &&
    	    (dis.atom2().isOxygen || dis.atom2().isNitrogen || dis.atom2().isSulfur)) {
    		// Calculating the HB sigmoid of atom1
    		Sigma.sigma(dis.distance(),
    					parameters.HBend[TsaiAtomicType1][TsaiAtomicType2],
    					parameters.HBp1[TsaiAtomicType1][TsaiAtomicType2],
    					parameters.HBp2[TsaiAtomicType1][TsaiAtomicType2],
    					parameters.HBvalAtp1[TsaiAtomicType1][TsaiAtomicType2],
    					parameters.HBvalAtp2[TsaiAtomicType1][TsaiAtomicType2]);
    		tmpSigmaHBatom1 = Sigma.s;
    		tmpSigma_TagHBatom1 = Sigma.s_tag;
    		// Calculating the HB sigmoid of atom2
    		Sigma.sigma(dis.distance(),
    					parameters.HBend[TsaiAtomicType2][TsaiAtomicType1],
    					parameters.HBp1[TsaiAtomicType2][TsaiAtomicType1],
    					parameters.HBp2[TsaiAtomicType2][TsaiAtomicType1],
    					parameters.HBvalAtp1[TsaiAtomicType2][TsaiAtomicType1],
    					parameters.HBvalAtp2[TsaiAtomicType2][TsaiAtomicType1]);
    		tmpSigmaHBatom2 = Sigma.s;
    		tmpSigma_TagHBatom2 = Sigma.s_tag;
    		// Calculating the PHB sigmoid of atom1
    		Sigma.sigma(dis.distance(),
    					parameters.HBend[TsaiAtomicType1][TsaiAtomicType2],
    					parameters.HBp1[TsaiAtomicType1][TsaiAtomicType2],
    					parameters.HBp2[TsaiAtomicType1][TsaiAtomicType2],
    					parameters.HBvalAtp1[TsaiAtomicType1][TsaiAtomicType2],
    					parameters.HBvalAtp2[TsaiAtomicType1][TsaiAtomicType2]);
    		tmpSigmaPHBatom1 = Sigma.s;
    		tmpSigma_TagPHBatom1 = Sigma.s_tag;
    		// Calculating the PHB sigmoid of atom2
    		Sigma.sigma(dis.distance(),
    					parameters.HBend[TsaiAtomicType2][TsaiAtomicType1],
    					parameters.HBp1[TsaiAtomicType2][TsaiAtomicType1],
    					parameters.HBp2[TsaiAtomicType2][TsaiAtomicType1],
    					parameters.HBvalAtp1[TsaiAtomicType2][TsaiAtomicType1],
    					parameters.HBvalAtp2[TsaiAtomicType2][TsaiAtomicType1]);
    		tmpSigmaPHBatom2 = Sigma.s;
    		tmpSigma_TagPHBatom2 = Sigma.s_tag;
/*
if (updateA && ((dis.atom1().residueNumber()==165) && dis.atom1().name().equals("OG1") ||
(dis.atom2().residueNumber()==165) && dis.atom2().name().equals("OG1"))) {
System.out.println("\n" + tmpSigmaHBatom1 + " " + tmpSigmaHBatom2 + " " +  tmpSigmaPHBatom1
+ " " +  tmpSigmaPHBatom2 + "\n");
}
		if (parameters.HBend[TsaiAtomicType1][TsaiAtomicType2]>0.1)
				nirCount++;
*/			
    		if ((tmpSigmaHBatom2>0) || (tmpSigmaHBatom1>0)) {	// Some sort of a HB interaction
    			// case 1:   base---O...O---base    or    base---O...N---base
    			if (!baseAtom[lut[dis.atom1().number()]].isHydrogen && 
    			    !baseAtom[lut[dis.atom2().number()]].isHydrogen) {
    			    solvateHBAngle.updateAndEvaluateAtoms1234(
   			    		baseAtom[lut[dis.atom1().number()]],
   			    		dis.atom1(),
   			    		dis.atom2(),
   			    		baseAtom[lut[dis.atom2().number()]]);
   			    	sigmaValues.sigmHBa1 = tmpSigmaHBatom1 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmHBa1dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom1*dis.dDistanceDx();
   			    	sigmaValues.dsigmHBa1dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom1*dis.dDistanceDy();
   			    	sigmaValues.dsigmHBa1dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom1*dis.dDistanceDz();
	    			sigmaValues.dsigmHBa1dx2 = -sigmaValues.dsigmHBa1dx1;
   					sigmaValues.dsigmHBa1dy2 = -sigmaValues.dsigmHBa1dy1;
   					sigmaValues.dsigmHBa1dz2 = -sigmaValues.dsigmHBa1dz1;
   			    	sigmaValues.sigmHBa2 = tmpSigmaHBatom2 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmHBa2dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom2*dis.dDistanceDx();
   			    	sigmaValues.dsigmHBa2dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom2*dis.dDistanceDy();
   			    	sigmaValues.dsigmHBa2dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom2*dis.dDistanceDz();
	    			sigmaValues.dsigmHBa2dx2 = -sigmaValues.dsigmHBa1dx1;
   					sigmaValues.dsigmHBa2dy2 = -sigmaValues.dsigmHBa1dy1;
   					sigmaValues.dsigmHBa2dz2 = -sigmaValues.dsigmHBa1dz1;
    				if (!solvateHBAngle.zeroDerivative()) {
    					sigmaValues.hbAngleDerivativeNonZero = true;
   						sigmaValues.dsigmHBa1dx1 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx2();
   						sigmaValues.dsigmHBa1dy1 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmHBa1dz1 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmHBa1dx2 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmHBa1dy2 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmHBa1dz2 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz3();
   						sigmaValues.dsigmHBa2dx1 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx2();
   						sigmaValues.dsigmHBa2dy1 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmHBa2dz1 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmHBa2dx2 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmHBa2dy2 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmHBa2dz2 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz3();
   						sigmaValues.dsigmHBa1dx3 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx1();
   						sigmaValues.dsigmHBa1dy3 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmHBa1dz3 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmHBa1dx4 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmHBa1dy4 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmHBa1dz4 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz4();
   						sigmaValues.dsigmHBa2dx3 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx1();
  						sigmaValues.dsigmHBa2dy3 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmHBa2dz3 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmHBa2dx4 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmHBa2dy4 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmHBa2dz4 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz4();
    				}
   			    	sigmaValues.sigmPHBa1 = tmpSigmaPHBatom1 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmPHBa1dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom1*dis.dDistanceDx();
   			    	sigmaValues.dsigmPHBa1dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom1*dis.dDistanceDy();
   			    	sigmaValues.dsigmPHBa1dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom1*dis.dDistanceDz();
	    			sigmaValues.dsigmPHBa1dx2 = -sigmaValues.dsigmPHBa1dx1;
   					sigmaValues.dsigmPHBa1dy2 = -sigmaValues.dsigmPHBa1dy1;
   					sigmaValues.dsigmPHBa1dz2 = -sigmaValues.dsigmPHBa1dz1;
   			    	sigmaValues.sigmPHBa2 = tmpSigmaPHBatom2 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmPHBa2dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom2*dis.dDistanceDx();
   			    	sigmaValues.dsigmPHBa2dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom2*dis.dDistanceDy();
   			    	sigmaValues.dsigmPHBa2dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom2*dis.dDistanceDz();
	    			sigmaValues.dsigmPHBa2dx2 = -sigmaValues.dsigmPHBa1dx1;
   					sigmaValues.dsigmPHBa2dy2 = -sigmaValues.dsigmPHBa1dy1;
   					sigmaValues.dsigmPHBa2dz2 = -sigmaValues.dsigmPHBa1dz1;
    				if (!solvateHBAngle.zeroDerivative()) {
   						sigmaValues.dsigmPHBa1dx1 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx2();
   						sigmaValues.dsigmPHBa1dy1 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmPHBa1dz1 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmPHBa1dx2 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmPHBa1dy2 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmPHBa1dz2 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz3();
   						sigmaValues.dsigmPHBa2dx1 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx2();
   						sigmaValues.dsigmPHBa2dy1 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmPHBa2dz1 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmPHBa2dx2 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmPHBa2dy2 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmPHBa2dz2 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz3();
   						sigmaValues.dsigmPHBa1dx3 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx1();
   						sigmaValues.dsigmPHBa1dy3 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmPHBa1dz3 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmPHBa1dx4 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmPHBa1dy4 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmPHBa1dz4 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz4();
   						sigmaValues.dsigmPHBa2dx3 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx1();
  						sigmaValues.dsigmPHBa2dy3 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmPHBa2dz3 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmPHBa2dx4 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmPHBa2dy4 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmPHBa2dz4 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz4();
    				}
   				}
   				else  // case 2:   N---H...O---base 
    			if (baseAtom[lut[dis.atom1().number()]].isHydrogen && 
    			    !baseAtom[lut[dis.atom2().number()]].isHydrogen) {
    			    solvateHBAngle.updateAndEvaluateAtoms1234(
   			    		dis.atom1(),
   			    		baseAtom[lut[dis.atom1().number()]],
   			    		dis.atom2(),
   			    		baseAtom[lut[dis.atom2().number()]]);   			    		
   			    	sigmaValues.sigmHBa1 = tmpSigmaHBatom1 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmHBa1dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom1*dis.dDistanceDx();
   			    	sigmaValues.dsigmHBa1dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom1*dis.dDistanceDy();
   			    	sigmaValues.dsigmHBa1dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom1*dis.dDistanceDz();
	    			sigmaValues.dsigmHBa1dx2 = -sigmaValues.dsigmHBa1dx1;
   					sigmaValues.dsigmHBa1dy2 = -sigmaValues.dsigmHBa1dy1;
   					sigmaValues.dsigmHBa1dz2 = -sigmaValues.dsigmHBa1dz1;
   			    	sigmaValues.sigmHBa2 = tmpSigmaHBatom2 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmHBa2dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom2*dis.dDistanceDx();
   			    	sigmaValues.dsigmHBa2dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom2*dis.dDistanceDy();
   			    	sigmaValues.dsigmHBa2dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom2*dis.dDistanceDz();
	    			sigmaValues.dsigmHBa2dx2 = -sigmaValues.dsigmHBa1dx1;
   					sigmaValues.dsigmHBa2dy2 = -sigmaValues.dsigmHBa1dy1;
   					sigmaValues.dsigmHBa2dz2 = -sigmaValues.dsigmHBa1dz1;
    				if (!solvateHBAngle.zeroDerivative()) {
    					sigmaValues.hbAngleDerivativeNonZero = true;
   						sigmaValues.dsigmHBa1dx1 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx1();
   						sigmaValues.dsigmHBa1dy1 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmHBa1dz1 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmHBa1dx2 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmHBa1dy2 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmHBa1dz2 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz3();
   						sigmaValues.dsigmHBa2dx1 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx1();
   						sigmaValues.dsigmHBa2dy1 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmHBa2dz1 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmHBa2dx2 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmHBa2dy2 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmHBa2dz2 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz3();
   						sigmaValues.dsigmHBa1dx3 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx2();
   						sigmaValues.dsigmHBa1dy3 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmHBa1dz3 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmHBa1dx4 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmHBa1dy4 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmHBa1dz4 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz4();
   						sigmaValues.dsigmHBa2dx3 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx2();
  						sigmaValues.dsigmHBa2dy3 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmHBa2dz3 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmHBa2dx4 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmHBa2dy4 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmHBa2dz4 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz4();
    				}
   			    	sigmaValues.sigmPHBa1 = tmpSigmaPHBatom1 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmPHBa1dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom1*dis.dDistanceDx();
   			    	sigmaValues.dsigmPHBa1dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom1*dis.dDistanceDy();
   			    	sigmaValues.dsigmPHBa1dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom1*dis.dDistanceDz();
	    			sigmaValues.dsigmPHBa1dx2 = -sigmaValues.dsigmPHBa1dx1;
   					sigmaValues.dsigmPHBa1dy2 = -sigmaValues.dsigmPHBa1dy1;
   					sigmaValues.dsigmPHBa1dz2 = -sigmaValues.dsigmPHBa1dz1;
   			    	sigmaValues.sigmPHBa2 = tmpSigmaPHBatom2 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmPHBa2dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom2*dis.dDistanceDx();
   			    	sigmaValues.dsigmPHBa2dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom2*dis.dDistanceDy();
   			    	sigmaValues.dsigmPHBa2dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom2*dis.dDistanceDz();
	    			sigmaValues.dsigmPHBa2dx2 = -sigmaValues.dsigmPHBa1dx1;
   					sigmaValues.dsigmPHBa2dy2 = -sigmaValues.dsigmPHBa1dy1;
   					sigmaValues.dsigmPHBa2dz2 = -sigmaValues.dsigmPHBa1dz1;
    				if (!solvateHBAngle.zeroDerivative()) {
   						sigmaValues.dsigmPHBa1dx1 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx1();
   						sigmaValues.dsigmPHBa1dy1 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmPHBa1dz1 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmPHBa1dx2 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmPHBa1dy2 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmPHBa1dz2 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz3();
   						sigmaValues.dsigmPHBa2dx1 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx1();
   						sigmaValues.dsigmPHBa2dy1 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmPHBa2dz1 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmPHBa2dx2 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmPHBa2dy2 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmPHBa2dz2 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz3();
   						sigmaValues.dsigmPHBa1dx3 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx2();
   						sigmaValues.dsigmPHBa1dy3 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmPHBa1dz3 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmPHBa1dx4 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmPHBa1dy4 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmPHBa1dz4 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz4();
   						sigmaValues.dsigmPHBa2dx3 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx2();
  						sigmaValues.dsigmPHBa2dy3 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmPHBa2dz3 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmPHBa2dx4 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmPHBa2dy4 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmPHBa2dz4 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz4();
    				}
   				}
   				else // case 3:   base---O...H---N 
    			if (!baseAtom[lut[dis.atom1().number()]].isHydrogen && 
    			    baseAtom[lut[dis.atom2().number()]].isHydrogen) {
    			    solvateHBAngle.updateAndEvaluateAtoms1234(
   			    		baseAtom[lut[dis.atom1().number()]],
   			    		dis.atom1(),
   			    		baseAtom[lut[dis.atom2().number()]],
   			    		dis.atom2());
   			    	sigmaValues.sigmHBa1 = tmpSigmaHBatom1 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmHBa1dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom1*dis.dDistanceDx();
   			    	sigmaValues.dsigmHBa1dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom1*dis.dDistanceDy();
   			    	sigmaValues.dsigmHBa1dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom1*dis.dDistanceDz();
	    			sigmaValues.dsigmHBa1dx2 = -sigmaValues.dsigmHBa1dx1;
   					sigmaValues.dsigmHBa1dy2 = -sigmaValues.dsigmHBa1dy1;
   					sigmaValues.dsigmHBa1dz2 = -sigmaValues.dsigmHBa1dz1;
   			    	sigmaValues.sigmHBa2 = tmpSigmaHBatom2 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmHBa2dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom2*dis.dDistanceDx();
   			    	sigmaValues.dsigmHBa2dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom2*dis.dDistanceDy();
   			    	sigmaValues.dsigmHBa2dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom2*dis.dDistanceDz();
	    			sigmaValues.dsigmHBa2dx2 = -sigmaValues.dsigmHBa1dx1;
   					sigmaValues.dsigmHBa2dy2 = -sigmaValues.dsigmHBa1dy1;
   					sigmaValues.dsigmHBa2dz2 = -sigmaValues.dsigmHBa1dz1;
    				if (!solvateHBAngle.zeroDerivative()) {
    					sigmaValues.hbAngleDerivativeNonZero = true;
   						sigmaValues.dsigmHBa1dx1 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx2();
   						sigmaValues.dsigmHBa1dy1 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmHBa1dz1 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmHBa1dx2 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmHBa1dy2 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmHBa1dz2 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz4();
   						sigmaValues.dsigmHBa2dx1 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx2();
   						sigmaValues.dsigmHBa2dy1 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmHBa2dz1 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmHBa2dx2 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmHBa2dy2 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmHBa2dz2 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz4();
   						sigmaValues.dsigmHBa1dx3 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx1();
   						sigmaValues.dsigmHBa1dy3 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmHBa1dz3 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmHBa1dx4 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmHBa1dy4 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmHBa1dz4 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz3();
   						sigmaValues.dsigmHBa2dx3 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx1();
  						sigmaValues.dsigmHBa2dy3 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmHBa2dz3 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmHBa2dx4 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmHBa2dy4 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmHBa2dz4 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz3();
    				}
   			    	sigmaValues.sigmPHBa1 = tmpSigmaPHBatom1 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmPHBa1dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom1*dis.dDistanceDx();
   			    	sigmaValues.dsigmPHBa1dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom1*dis.dDistanceDy();
   			    	sigmaValues.dsigmPHBa1dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom1*dis.dDistanceDz();
	    			sigmaValues.dsigmPHBa1dx2 = -sigmaValues.dsigmPHBa1dx1;
   					sigmaValues.dsigmPHBa1dy2 = -sigmaValues.dsigmPHBa1dy1;
   					sigmaValues.dsigmPHBa1dz2 = -sigmaValues.dsigmPHBa1dz1;
   			    	sigmaValues.sigmPHBa2 = tmpSigmaPHBatom2 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmPHBa2dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom2*dis.dDistanceDx();
   			    	sigmaValues.dsigmPHBa2dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom2*dis.dDistanceDy();
   			    	sigmaValues.dsigmPHBa2dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom2*dis.dDistanceDz();
	    			sigmaValues.dsigmPHBa2dx2 = -sigmaValues.dsigmPHBa1dx1;
   					sigmaValues.dsigmPHBa2dy2 = -sigmaValues.dsigmPHBa1dy1;
   					sigmaValues.dsigmPHBa2dz2 = -sigmaValues.dsigmPHBa1dz1;
    				if (!solvateHBAngle.zeroDerivative()) {
   						sigmaValues.dsigmPHBa1dx1 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx2();
   						sigmaValues.dsigmPHBa1dy1 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmPHBa1dz1 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmPHBa1dx2 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmPHBa1dy2 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmPHBa1dz2 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz4();
   						sigmaValues.dsigmPHBa2dx1 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx2();
   						sigmaValues.dsigmPHBa2dy1 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmPHBa2dz1 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmPHBa2dx2 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmPHBa2dy2 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmPHBa2dz2 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz4();
   						sigmaValues.dsigmPHBa1dx3 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx1();
   						sigmaValues.dsigmPHBa1dy3 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmPHBa1dz3 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmPHBa1dx4 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmPHBa1dy4 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmPHBa1dz4 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz3();
   						sigmaValues.dsigmPHBa2dx3 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx1();
  						sigmaValues.dsigmPHBa2dy3 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmPHBa2dz3 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmPHBa2dx4 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmPHBa2dy4 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmPHBa2dz4 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz3();
    				}
   				}
   				else // case 4:   N---H...H---N  
    			if (baseAtom[lut[dis.atom1().number()]].isHydrogen && 
    			    baseAtom[lut[dis.atom2().number()]].isHydrogen) {
    			    solvateHBAngle.updateAndEvaluateAtoms1234(
   			    		dis.atom1(),
   			    		baseAtom[lut[dis.atom1().number()]],
   			    		baseAtom[lut[dis.atom2().number()]],
   			    		dis.atom2());
   			    	sigmaValues.sigmHBa1 = tmpSigmaHBatom1 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmHBa1dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom1*dis.dDistanceDx();
   			    	sigmaValues.dsigmHBa1dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom1*dis.dDistanceDy();
   			    	sigmaValues.dsigmHBa1dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom1*dis.dDistanceDz();
	    			sigmaValues.dsigmHBa1dx2 = -sigmaValues.dsigmHBa1dx1;
   					sigmaValues.dsigmHBa1dy2 = -sigmaValues.dsigmHBa1dy1;
   					sigmaValues.dsigmHBa1dz2 = -sigmaValues.dsigmHBa1dz1;
   			    	sigmaValues.sigmHBa2 = tmpSigmaHBatom2 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmHBa2dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom2*dis.dDistanceDx();
   			    	sigmaValues.dsigmHBa2dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom2*dis.dDistanceDy();
   			    	sigmaValues.dsigmHBa2dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagHBatom2*dis.dDistanceDz();
	    			sigmaValues.dsigmHBa2dx2 = -sigmaValues.dsigmHBa1dx1;
   					sigmaValues.dsigmHBa2dy2 = -sigmaValues.dsigmHBa1dy1;
   					sigmaValues.dsigmHBa2dz2 = -sigmaValues.dsigmHBa1dz1;
    				if (!solvateHBAngle.zeroDerivative()) {
    					sigmaValues.hbAngleDerivativeNonZero = true;
   						sigmaValues.dsigmHBa1dx1 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx1();
   						sigmaValues.dsigmHBa1dy1 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmHBa1dz1 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmHBa1dx2 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmHBa1dy2 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmHBa1dz2 += tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz4();
   						sigmaValues.dsigmHBa2dx1 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx1();
   						sigmaValues.dsigmHBa2dy1 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmHBa2dz1 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmHBa2dx2 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmHBa2dy2 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmHBa2dz2 += tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz4();
   						sigmaValues.dsigmHBa1dx3 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx2();
   						sigmaValues.dsigmHBa1dy3 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmHBa1dz3 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmHBa1dx4 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmHBa1dy4 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmHBa1dz4 = tmpSigmaHBatom1*solvateHBAngle.DhbAngScoreDz3();
   						sigmaValues.dsigmHBa2dx3 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx2();
  						sigmaValues.dsigmHBa2dy3 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmHBa2dz3 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmHBa2dx4 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmHBa2dy4 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmHBa2dz4 = tmpSigmaHBatom2*solvateHBAngle.DhbAngScoreDz3();
    				}
   			    	sigmaValues.sigmPHBa1 = tmpSigmaPHBatom1 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmPHBa1dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom1*dis.dDistanceDx();
   			    	sigmaValues.dsigmPHBa1dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom1*dis.dDistanceDy();
   			    	sigmaValues.dsigmPHBa1dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom1*dis.dDistanceDz();
	    			sigmaValues.dsigmPHBa1dx2 = -sigmaValues.dsigmPHBa1dx1;
   					sigmaValues.dsigmPHBa1dy2 = -sigmaValues.dsigmPHBa1dy1;
   					sigmaValues.dsigmPHBa1dz2 = -sigmaValues.dsigmPHBa1dz1;
   			    	sigmaValues.sigmPHBa2 = tmpSigmaPHBatom2 * solvateHBAngle.hbAngScore();
   			    	sigmaValues.dsigmPHBa2dx1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom2*dis.dDistanceDx();
   			    	sigmaValues.dsigmPHBa2dy1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom2*dis.dDistanceDy();
   			    	sigmaValues.dsigmPHBa2dz1 = solvateHBAngle.hbAngScore()*tmpSigma_TagPHBatom2*dis.dDistanceDz();
	    			sigmaValues.dsigmPHBa2dx2 = -sigmaValues.dsigmPHBa1dx1;
   					sigmaValues.dsigmPHBa2dy2 = -sigmaValues.dsigmPHBa1dy1;
   					sigmaValues.dsigmPHBa2dz2 = -sigmaValues.dsigmPHBa1dz1;
    				if (!solvateHBAngle.zeroDerivative()) {
   						sigmaValues.dsigmPHBa1dx1 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx1();
   						sigmaValues.dsigmPHBa1dy1 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmPHBa1dz1 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmPHBa1dx2 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmPHBa1dy2 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmPHBa1dz2 += tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz4();
   						sigmaValues.dsigmPHBa2dx1 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx1();
   						sigmaValues.dsigmPHBa2dy1 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy1();
   						sigmaValues.dsigmPHBa2dz1 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz1();
   						sigmaValues.dsigmPHBa2dx2 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx4();
   						sigmaValues.dsigmPHBa2dy2 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy4();
   						sigmaValues.dsigmPHBa2dz2 += tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz4();
   						sigmaValues.dsigmPHBa1dx3 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx2();
   						sigmaValues.dsigmPHBa1dy3 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmPHBa1dz3 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmPHBa1dx4 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmPHBa1dy4 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmPHBa1dz4 = tmpSigmaPHBatom1*solvateHBAngle.DhbAngScoreDz3();
   						sigmaValues.dsigmPHBa2dx3 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx2();
  						sigmaValues.dsigmPHBa2dy3 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy2();
   						sigmaValues.dsigmPHBa2dz3 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz2();
   						sigmaValues.dsigmPHBa2dx4 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDx3();
   						sigmaValues.dsigmPHBa2dy4 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDy3();
   						sigmaValues.dsigmPHBa2dz4 = tmpSigmaPHBatom2*solvateHBAngle.DhbAngScoreDz3();
    				}
   				}
   			} // Of checkin if a HB interaction soes exist
   		} // Of checking for HB forming pair 		
   	} // of updateSigmVals



/** 
 * For any polar atom (O,N or S in Cys) we calculate the base atom that participate in the definition of the hydrogen
 * bond angle. This base atom is the attached hydrogen (if present), or the heavy atom to which the polar atom is 
 * attached (when the hydrogen is not present).
 **/
	private final void setBaseAtom() {
		Atom atom1,atom2;
		for (int c1=0 ; c1<atomList.size() ; c1++) {
			atom1 = atomList.atomAt(c1);
			for (int c2=0 ; c2<atom1.bonded().size() ; c2++) {
				atom2 = atom1.bonded().atomAt(c2);
				// Treating OXYGEN atoms. 
				// This is easy because the oxygens in proteins 
				// are always tied to one atom only.
				if (atom1.isOxygen)
					baseAtom[lut[atom1.number()]] = atom2;
				// Treating NITROGEN atoms.
				// First, the case of amides in glutamines and asparagines
				if (!atom2.isHydrogen && ((atom1.type==NND) || (atom1.type==QNE)))
					baseAtom[lut[atom1.number()]] = atom2;
				// Second, the case of amides without explicit H attached
				if ((atom1.type==KNZ) || (atom1.type==RNH) || (atom1.type==TRN))
					baseAtom[lut[atom1.number()]] = atom2;
				// Third , regular H attached
				if (atom2.isHydrogen && ((atom1.type==HND) || (atom1.type==HNE) || 
				    (atom1.type==RNE) || (atom1.type==WNE) ||
					(atom1.type==AN) ||
					(atom1.type==CN) ||
					(atom1.type==DN) ||
					(atom1.type==EN) ||
					(atom1.type==FN) ||
					(atom1.type==GN) ||
					(atom1.type==HN) ||
					(atom1.type==IN) ||
					(atom1.type==KN) ||
					(atom1.type==LN) ||
					(atom1.type==MN) ||
					(atom1.type==NN) ||
					(atom1.type==QN) ||
					(atom1.type==RN) ||
					(atom1.type==SN) ||
					(atom1.type==TN) ||
					(atom1.type==VN) ||
					(atom1.type==WN) ||
					(atom1.type==YN)))
					baseAtom[lut[atom1.number()]] = atom2;
				// Treating the SULFUR atoms of Cystines.
				if (atom1.type==CSG)
					baseAtom[lut[atom1.number()]] = atom2;
			}
		}
	}
    
}
	
