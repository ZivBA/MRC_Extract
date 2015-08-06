/*
 * Created on 17/11/2004
 * Window - Preferences - Java - Code Style - Code Templates
 */
package meshi.energy.hydrogenBond;

import java.util.Iterator;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.Protein;
import meshi.parameters.AtomTypes;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.MeshiException;

/**
 * @author amilev
 *
 * This class is used to create a HydroenBondEnergy.
 * The user does not need to know how to create such an energy term but only to give its weight. 
 */
public class HydrogenBondsCreator extends EnergyCreator implements KeyWords, AtomTypes {

    //------------------------- data fields ---------------------------
    
    //private IsN isN = new IsN(); old version
    HydrogenBondsEnergy hydrogenBondsEnergy;


    /**
  * @return Returns the hydrogenBondsEnergy.
  */
 public final HydrogenBondsEnergy getHydrogenBondsEnergy() {	return hydrogenBondsEnergy;	}

    private  DistanceList specialDis = null;
    public final DistanceList getSpecialDis(){return specialDis ;}

    int[ ] specialDisArray = null;
    public int[] getSpecialDisArray() {
        return specialDisArray;
    }

    private boolean antiParalel;
    public boolean getAntiParalel() {
        return antiParalel;
    }




    //------------------------- constructors ---------------------------

    public HydrogenBondsCreator() {
		super(HYDROGEN_BONDS);
	}
    
	/**
	 * @param weight
	 */
	public HydrogenBondsCreator(double weight) {
		super(weight);
	}

	public HydrogenBondsCreator(double weight , int[ ] specialDisArray,Protein protein,DistanceMatrix distanceMatrix,boolean antiParallel ){
		super(weight);
        readSpecialDistance(specialDisArray ,protein,distanceMatrix,antiParallel);
    }

     public HydrogenBondsCreator(int[ ] specialDisArray,Protein protein,DistanceMatrix distanceMatrix,boolean antiParallel){
		super(HYDROGEN_BONDS);
        readSpecialDistance(specialDisArray ,protein,distanceMatrix,antiParallel);
    }

    private void readSpecialDistance(int[ ] specialDisArray,Protein protein,DistanceMatrix distanceMatrix, boolean antiParallel){
        this.antiParalel = antiParallel;
        if(specialDisArray .length % 2 != 0)
                        throw new MeshiException("HydrogenBondCreator Should get array of even length");
        this.specialDisArray = specialDisArray ;
        int length = specialDisArray.length;
             specialDis = new DistanceList();
        if(antiParallel){
            for (int i = 0; i<length-1; i=i+2){
                specialDis .fastAdd(distanceMatrix.distance(protein.residue(specialDisArray [i]).getAtom("H"),protein.residue(specialDisArray [i+1]).getAtom("O") ) );
                specialDis .fastAdd(distanceMatrix.distance(protein.residue(specialDisArray [i+1]).getAtom("H"),protein.residue(specialDisArray [i]).getAtom("O") ) );
            }
        }
        else{
            for(int i=0;i<length-1;i=i+2){
                specialDis .fastAdd(distanceMatrix.distance(protein.residue(specialDisArray [i]).getAtom("O"), protein.residue(specialDisArray [i+1]+1).getAtom("H")));
                specialDis .fastAdd(distanceMatrix.distance(protein.residue(specialDisArray [i]).getAtom("H"), protein.residue(specialDisArray [i+1]- 1).getAtom("O") ) );
            }
        }
            specialDis .print();
    }



    //------------------------------ methods ------------------------------

    /*
     * @param protein does not been used; must be given to implement an abstruct method of EnergyCreator.
     * @param distanceMatrix
     * @param commands gives the path to the data directory
     */
	public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
                                           CommandList commands) {
        //create parameters
		if (parametersList== null)
            {
                //AtomList nitrogens = protein.atoms().filter(isN); old version
                parametersList = new HydrogenBondsParametersList(parametersDirectory(commands)+
                                                                 "/"+LENNARD_JONES_PARAMETERS);//,nitrogens);
            }
        Iterator atoms = distanceMatrix.atomList().iterator();
        Atom atom;
        while ((atom = (Atom) atoms.next()) != null) {
            if (IsHO.isH(atom)) atom.addAttribute(new HB_AtomAttribute(true,false));
            if (IsHO.isO(atom)) atom.addAttribute(new HB_AtomAttribute(false,true));
        }

        if(specialDis != null ){

           hydrogenBondsEnergy = new HydrogenBondsEnergy(distanceMatrix,
                                                                  (HydrogenBondsParametersList) parametersList,
                                                                  weight(),
                                                                  new HBondList(distanceMatrix,(HydrogenBondsParametersList)parametersList,specialDis ),
                                                                  specialDis);

        }
        else{
        hydrogenBondsEnergy = new HydrogenBondsEnergy(distanceMatrix,
                                                      ( HydrogenBondsParametersList) parametersList,
                                                      weight(),
                                                      new HBondList(distanceMatrix,(HydrogenBondsParametersList)parametersList));
        }
        distanceMatrix.energyTermsDistanceLists() .fastAdd(HBondList .inputNewHBList() );
        return hydrogenBondsEnergy;
	}
	
	
    //---------------------------------------------------------------------------
	/* private static class IsN implements Filter,AtomTypes{
       int[] sortBB_NITROGENS;
        
       //create a new SORT copy of BB_NITROGENS that are defined in interface AtomType
       public IsN(){             
       sortBB_NITROGENS = new int[BB_NITROGENS.length];
       for(int i = 0;i<BB_NITROGENS.length;i++)
       sortBB_NITROGENS[i] = BB_NITROGENS[i];
       Arrays.sort(sortBB_NITROGENS);
       }
        
       public boolean accept(Object obj) {
       return (Arrays.binarySearch(sortBB_NITROGENS,((Atom)obj).type) >= 0);
            
       }
       }//isN
    */
}
