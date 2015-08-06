package meshi.molecularElements.hydrogenBonds;

import meshi.energy.solvate.SolvateParametersList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.parameters.AtomTypes;

/** 
 * This class describe a list of hydrogen bonds in a protein. The hydrogen bonds are described as in:
 * Dahiyat et al. Protein Sci. 1996 May;5(5):895-903. See the "DahiyatHydrogenBond" class for a much more
 * detailed description.
 * 
 * The construcor requires a parameter that points to an instance of "SolvateParametersList". From this instance 
 * the distance dependence of the hydrogen bonds are extracted. See the heading documentation of the 
 * SolvateParametersList class for more details. 
 * 
 * 
 * This parameter 30 determines the refresh rate of the HB vector. Once every this number of updates, broken hydrogen 
 * bonds that are no longer in the non-bonded-list are removed from the vector. This number should be >>1 so that 
 * the refresh does not impend the updates too much, but 30 was just a thumb figure.   	
 *
 **/

public class HydrogenBondDahiyatList extends AbstractHydrogenBondList implements AtomTypes {
	
	/**
	 * This array stores pointers to the base atoms of every non-hydrogen polar atom in the protein. 
	 * The indexation is through the atom number (field)   
	 **/
    protected Atom[] baseAtom; 
	/**
	 * Since currently the hydrogen bond list is only used by the solvate energy, we use the parameters 
	 * for the hydrogen bond from this class. Later, this can be changed.   
	 **/
    protected SolvateParametersList parameters;
    
    /**
     * In this class are the angle parameters for the HB
     **/
    protected DahiyatImplementationConstants angleParameters = new DahiyatHighAccuracyAngleParamaters();


    public HydrogenBondDahiyatList() {
    	throw new RuntimeException("\nERROR: without parameters the hydrogen bonds cannot be formed.\n");    	
    }

    public HydrogenBondDahiyatList(DistanceMatrix dm, AtomList atomList,SolvateParametersList parameters) {
    	super(dm, atomList, 30 /* DEFAULT_REFRESH_VECTOR_RATE */);
    	this.parameters = parameters;
    	try {
    		update(false,0);
    	}
    	catch (Exception e) {
    		System.out.print("\nAn error occur while creating the Dahiyat hydrogen bond list.\n\n");
    		e.printStackTrace();
    		System.out.print("\n\n");
    		throw new RuntimeException("");
    	}
    }

    
		
    /**
      * This method updates the pointers in the baseAtom array. For non-polar atoms they should remain null.
      * For any polar atom (O,N or S in Cys) we calculate the base atom that participate in the definition of 
      * the hydrogen bond angle. This base atom is the attached hydrogen (if present), or the heavy atom to which 
      * the polar atom is attached (when the hydrogen is not present).      
      * 
      **/
    protected void buildSpecificStructures() {
		baseAtom = new Atom[maxAtomNum+1];
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

	/**
	 * Creating a Dahiyat-like hydrogen-bond.
	 **/
    protected AbstractHydrogenBond createHBfromPolars(Atom atom1,Atom atom2) {
		if ((atom1.type==MSD) || (atom2.type==MSD) || (atom1.type==PN) || (atom2.type==PN))
			return null;

    	AtomList tmpList = new AtomList();
    	tmpList.add(atom1);
    	tmpList.add(atom2);
    	tmpList.add(baseAtom[lut[atom1.number()]]);
    	tmpList.add(baseAtom[lut[atom2.number()]]);
    	
    	int TsaiAtomicType1 = parameters.atomicTypeConverter[atom1.type]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
    	int TsaiAtomicType2 = parameters.atomicTypeConverter[atom2.type]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
    	double end = parameters.HBend[TsaiAtomicType1][TsaiAtomicType2];
    	double p1 = parameters.HBp1[TsaiAtomicType1][TsaiAtomicType2];
    	double p2 = parameters.HBp2[TsaiAtomicType1][TsaiAtomicType2];
		double valAtp1 = parameters.HBvalAtp1[TsaiAtomicType1][TsaiAtomicType2];
    	double valAtp2 = parameters.HBvalAtp2[TsaiAtomicType1][TsaiAtomicType2];
    	
    	if (end<0.1)  // With the current paramters for the solvate energy: end=0 means that hydrogen bonding is not relevent to this pair.
			return null;
    	
    	return new HydrogenBondDahiyat(dm, tmpList,	 p1,  p2,  end,  valAtp1,  valAtp2, angleParameters );

    }


}


  
