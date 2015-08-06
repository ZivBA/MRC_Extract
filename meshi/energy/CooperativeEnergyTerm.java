package meshi.energy;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.UpdateableException;
/** 
 * A super class for those energy terms that operate on the entire protein, or on a specific group of atoms that
 * occur only once. 
 * Typical cases are Radius of Gyration, or solvation energies. 
 **/
public abstract class CooperativeEnergyTerm extends AbstractEnergy  {

    protected DistanceMatrix dm;
    /** 
     *The atoms relevent to this energy.
     **/
    protected AtomList atomList; 

    public CooperativeEnergyTerm() {}

    public CooperativeEnergyTerm(Object[] updateableResources,
                                 AtomList atomList, 
                                 DistanceMatrix dm,
                                 ParametersList parameters, 
                                 double weight) {
	super(addDistanceMatrixToUpdateable(updateableResources,dm) , parameters , weight);
	this.dm = dm;
	this.atomList = atomList;
    }


    /**
     * Test the accuracy of deriving the energy function by coordinate number "i".
     * The analytical derivation is compared to the numerical one.
     **/
    public void test(TotalEnergy totalEnergy,Atom atom){
        if (atomList == null) 
	    throw new RuntimeException("Cannot test "+this+"\n"+"No atoms defined");
        if(atomList.whereIs(atom) < 0)
            return;

        double[][] coordinates = new double[3][];
        coordinates[0] = atom.X();
        coordinates[1] = atom.Y();
        coordinates[2] = atom.Z();
        for(int i = 0; i< 3; i++) {
            try{totalEnergy.update();}catch(UpdateableException ue){}
            double x = coordinates[i][0];
            coordinates[i][1] = 0;
            double e1 = evaluate();
            double analiticalForce = coordinates[i][1];
            coordinates[i][0] += DX;
            // Whatever should be updated ( such as distance matrix torsion list etc. )
            try{totalEnergy.update();}catch(UpdateableException ue){}
            double e2 = evaluate();
            double de = e2-e1;
            double numericalForce = - de/DX;
            coordinates[i][0] -= DX;
            double diff = Math.abs(analiticalForce - numericalForce);
            
            if ((2*diff/(Math.abs(analiticalForce)+Math.abs(numericalForce)+verySmall)) > relativeDiffTolerance){
                System.out.println("Testing "+this);
                System.out.println("Atom["+atom.number()+"]."+XYZ.charAt(i)+" = "+x);
                System.out.println("Analytical force = "+analiticalForce);
                System.out.println("Numerical force  = "+numericalForce);
                
                System.out.println("diff = "+diff+"\n"+
                                   "tolerance = 2*diff/(|analiticalForce| + |numericalForce|+verySmall) = "+
                                   2*diff/(Math.abs(analiticalForce) + Math.abs(numericalForce)+verySmall));
                System.out.println();
            }
            if ((e1 == AbstractEnergy.INFINITY) | (e1 == AbstractEnergy.NaN))
                System.out.println("Testing "+this+"\ne1 = "+e1);
            if ((e2 == AbstractEnergy.INFINITY) | (e2 == AbstractEnergy.NaN))
                System.out.println("Testing "+this+"\ne2 = "+e2);
            if ((analiticalForce == AbstractEnergy.INFINITY) | (analiticalForce == AbstractEnergy.NaN))
                System.out.println("Testing "+this+"\nanaliticalForce = "+analiticalForce);

        }
    }

    /**
     * This auxillary method will add the DistanceMatrix object, dm, to the object list, unless there is already an instance of DistanceMatrix there.
     * It also put the DistanceMatrix object first.
     **/
    protected static Object[] addDistanceMatrixToUpdateable(Object[] ar, DistanceMatrix dm) {
    	if (ar == null)
    		return toArray(dm);
    	boolean dmInArray = false;
    	int indOfDM = -1;
    	for (int c=0 ; c<ar.length ; c++)          
    		if (ar[c] instanceof DistanceMatrix) { 
    			dmInArray = true;
    			indOfDM = c;
    		}
    	Object[] newAr;
    	if (dmInArray) {
    		newAr = new Object[ar.length];
    		for (int c=0 ; c<ar.length ; c++) 
    			newAr[c] = ar[c];
    		Object tmp = newAr[0];
    		newAr[0] = newAr[indOfDM];
    		newAr[indOfDM] = tmp;
    	}
    	else {
    		newAr = new Object[ar.length+1];
    		newAr[0] = dm;
    		for (int c=0 ; c<ar.length ; c++) 
    			newAr[c+1] = ar[c];
    	}
    	return newAr;
    }

}


  
