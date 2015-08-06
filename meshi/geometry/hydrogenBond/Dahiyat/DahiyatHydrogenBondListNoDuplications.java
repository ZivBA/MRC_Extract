package meshi.geometry.hydrogenBond.Dahiyat;

import meshi.geometry.DistanceMatrix;
import meshi.geometry.hydrogenBond.AbstractHydrogenBond;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.parameters.Residues;
import meshi.util.UpdateableException;
import meshi.util.mathTools.HB_Sigma;

/**
 * 
 * The purpose of this class is to make sure that with a given polar atom only one atom from a ASP GLU or ARG sidechain will make 
 * a hydrogen bond. If another hydrogen bond is formed with the other atom in the functional group the weaker of the two will have
 * it's value set to zero. 
 * 
 * @author Nir
 *
 */

public class DahiyatHydrogenBondListNoDuplications extends
		DahiyatHydrogenBondList implements Residues {

    protected Atom[] lutSameGroup; // If i is the atom number of an atom in the funcional group of ASP GLU or ARG, then lutSameGroup[i] is a pointer to the other atom in the functional group. It is null for all other atoms.     

	
	public DahiyatHydrogenBondListNoDuplications() {
    	throw new RuntimeException("\nERROR: without parameters the hydrogen bonds cannot be formed.\n");    	
	}

	public DahiyatHydrogenBondListNoDuplications(DistanceMatrix dm,
			AtomList atomList, DahiyatParametersInterface parameters) {
		super(dm, atomList, parameters);
	}

    protected DahiyatHydrogenBond createSpecificBond(AtomList tmpList,	HB_Sigma tmpSigma) {
    	return new DahiyatHydrogenBondNoDuplications(dm, tmpList,	tmpSigma, parameters);
    }
	
    /**
     * In addition to updating the regular parmeters to the Dahiyat list, this method also find
     * the ARG, ASP and GLU groups, where duplications of hydrogen bonds can occur.
     */
    protected void buildSpecificStructures() {
    	super.buildSpecificStructures();
    	lutSameGroup = new Atom[maxAtomNum+1];
		Atom otherAtom = null;
		boolean functionalGroup = false;
    	for (int c=0; c<atomList.size() ; c++) {
    		functionalGroup = false;
    		otherAtom = null;
    		if ((atomList.atomAt(c).residue().type==ASP) && atomList.atomAt(c).name().equals("OD1")) {
    			otherAtom = atomList.findAtomInList("OD2", atomList.atomAt(c).residueNumber());
    			functionalGroup = true;
    		}
    		if ((atomList.atomAt(c).residue().type==ASP) && atomList.atomAt(c).name().equals("OD2")) {
    			otherAtom = atomList.findAtomInList("OD1", atomList.atomAt(c).residueNumber());
    			functionalGroup = true;
    		}
    		if ((atomList.atomAt(c).residue().type==GLU) && atomList.atomAt(c).name().equals("OE1")) {
    			otherAtom = atomList.findAtomInList("OE2", atomList.atomAt(c).residueNumber());
    			functionalGroup = true;
    		}
    		if ((atomList.atomAt(c).residue().type==GLU) && atomList.atomAt(c).name().equals("OE2")) {
    			otherAtom = atomList.findAtomInList("OE1", atomList.atomAt(c).residueNumber());
    			functionalGroup = true;
    		}
    		if ((atomList.atomAt(c).residue().type==ARG) && atomList.atomAt(c).name().equals("NH1")) {
    			otherAtom = atomList.findAtomInList("NH2", atomList.atomAt(c).residueNumber());
    			functionalGroup = true;
    		}
    		if ((atomList.atomAt(c).residue().type==ARG) && atomList.atomAt(c).name().equals("NH2")) {
    			otherAtom = atomList.findAtomInList("NH1", atomList.atomAt(c).residueNumber());
    			functionalGroup = true;
    		}
    		if (functionalGroup) {
    			if (otherAtom!=null) 
    				lutSameGroup[atomList.atomAt(c).number()] = otherAtom;
    			else {
    				lutSameGroup[atomList.atomAt(c).number()] = atomList.atomAt(c); // Yes, in case of a stupid error like having one functional atom but not the other, we will make the atom its functional twin.
    				System.out.println("WARNING: The functional group in residue: " + atomList.atomAt(c).residue() + 
    						" consists of only one polar atom, where there should be 2. The run can continue, but this is a sign of a very poor model. Was the crystalographer drunk?");
    			}
    		}
    	}    	
    }
    
    
    public void update(int updateNumber) throws UpdateableException {
    	super.update(updateNumber);
    	// Now updating the duplication weights in each HB
    	DahiyatHydrogenBondNoDuplications hb;
    	DahiyatHydrogenBondNoDuplications otherHB;
    	double value1, value2;
        for(AbstractHydrogenBond nonCastedHB : bondList)
        	if (nonCastedHB.active()) {
        		hb = (DahiyatHydrogenBondNoDuplications) nonCastedHB;
        		if (hb.isPolar1Relevant()) {
        			value1 = hb.hbValWithoutDuplicationWeight();
        			otherHB = (DahiyatHydrogenBondNoDuplications) findBondByPolars(lutSameGroup[hb.getFirstPolar().number()], hb.getSecondPolar());
        			if (otherHB!=null) {
        				value2 = otherHB.hbValWithoutDuplicationWeight();
        				if (value2>value1) {
        					hb.setDuplicationWeight(0.0);
        					otherHB.setDuplicationWeight(1.0);
        					//System.out.println(hb + "  " + value1 + "\n" +  otherHB + "  " + value2 + "\n\n");
        				}
        				else {
        					otherHB.setDuplicationWeight(0.0);
        					hb.setDuplicationWeight(1.0);
        					//System.out.println(hb + "  " + value1 + "\n" +  otherHB + "  " + value2 + "\n\n");
        				}
        			}
        		}
        		if (hb.isPolar2Relevant()) {
        			value1 = hb.hbValWithoutDuplicationWeight();
        			otherHB = (DahiyatHydrogenBondNoDuplications) findBondByPolars(lutSameGroup[hb.getSecondPolar().number()], hb.getFirstPolar());
        			if (otherHB!=null) {
        				value2 = otherHB.hbValWithoutDuplicationWeight();
        				if (value2>value1) {
        					hb.setDuplicationWeight(0.0);
        					otherHB.setDuplicationWeight(1.0);
        					//System.out.println(hb + "  " + value1 + "\n" +  otherHB + "  " + value2 + "\n\n");
        				}
        				else {
        					otherHB.setDuplicationWeight(0.0);
        					hb.setDuplicationWeight(1.0);
        					//System.out.println(hb + "  " + value1 + "\n" +  otherHB + "  " + value2 + "\n\n");
        				}
        			}
        		}
        	}
    }
    
}
