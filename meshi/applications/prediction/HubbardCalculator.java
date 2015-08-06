package meshi.applications.prediction;
import meshi.applications.hubbard.HubbardPosition;
import meshi.molecularElements.AtomList;

public class HubbardCalculator 
{
	public static double[][] C1;
	public static double[][] C2;
	public static HubbardPosition[] Popo;
	public static int len;
    public static final int sub = 4; //the length of the first sub protein that we check.
	                   //According to Hubbard plot - suppose to be 3 but there is a problem 
	                   //to the Hubbard overlap to handle it.
    public static final double maximal_distance = 4.9;//the maximal distance between 2 Ca atoms for joining them to the
                                       //group of Atoms we send at the next iteration.
    public static final int num_of_loops = 20;//the maximal number of iterations over the protein in finding
                                    //the best superposition.


    public static double hubbard_gdt(AtomList reference, AtomList protein1) {
    	return hubbard_gdt(reference, protein1, 1, 2, 4, 8);
    }
    
    public static double hubbard_gdt(AtomList reference, AtomList protein1 , double cutoff1, double cutoff2, double cutoff3, double cutoff4) {
        //creating two arrays (one for each given structure)
        //according to the atoms coordinates writen at the files.
        read_files(reference, protein1);
        Popo = new HubbardPosition[len+1-sub];//represents all the superspositions of the protein.
        initialize(Popo);//initialize Popo (for each subunit of the 
                                                         //protein at the base subunit size.
        loop_over_Popo(Popo);//for each start point (according to the base subunit)
                                                    //finding the best superposition of the protein.
        double[] best_distance = new double[len];
        build_best_distance(best_distance, Popo);//building the best distances array.
System.out.println(cutoff1 + " " + fracBelowRMS(best_distance,cutoff1,reference.CAFilter().size()));
System.out.println(cutoff2 + " " + fracBelowRMS(best_distance,cutoff2,reference.CAFilter().size()));
System.out.println(cutoff3 + " " + fracBelowRMS(best_distance,cutoff3,reference.CAFilter().size()));
System.out.println(cutoff4 + " " + fracBelowRMS(best_distance,cutoff4,reference.CAFilter().size()));
        return (0.25*fracBelowRMS(best_distance,cutoff1,reference.CAFilter().size()) +
        		0.25*fracBelowRMS(best_distance,cutoff2,reference.CAFilter().size()) +
        		0.25*fracBelowRMS(best_distance,cutoff3,reference.CAFilter().size()) +
        		0.25*fracBelowRMS(best_distance,cutoff4,reference.CAFilter().size()));
       
    }//gdt

 
    //---------------------------------------------------------------------------------------
    
    private static void initialize(HubbardPosition[] Popo) {
    	int place;
    	for(place = 0; place< Popo.length; place++)   //initialize the basic states
            {
                Popo[place] = new HubbardPosition(sub, place, maximal_distance, C1 , C2);
            }
    }//End of initializing

     
     //--------------------------------------------------------------------------------------
    private static void loop_over_Popo(HubbardPosition[] Popo) {
        for(int i = 0; i<Popo.length; i++)//go over all our HubbardPositions.
            {
                Popo[i].find_best_conformation(num_of_loops);	
            }//each position now has an array of the smallest sum of distances
    }

    
    //---------------------------------------------------------------------------------------
    
    private static void build_best_distance(double[] best_distance, HubbardPosition[] Popo) {
    	int i,j;
        for (i=0; i<len; i++) { // building the best_distance array
           	double temp_best = 100000000;
           	for(j=0; j<Popo.length; j++) {
           		if(Popo[j].numbers[i] < temp_best) {
           			temp_best = Popo[j].numbers[i];
           		}
           	}
           	best_distance[i] = Math.sqrt(temp_best);
        } //the best distance's array is ready
    }//build_best_distance

    
    //-----------------------------------------------------------------------------------------------
	// find the fraction of the protein overlap (length of len) that is below a certain rms
	// on the Hubbard plot curve. 
    private static double fracBelowRMS(double[] best_distance , double rms , int referLength) {
    	for (int c=0 ; c<len ; c++) 
    		if (best_distance[c]>rms) {
    			return c*1.0/referLength;
		}
    	return len*1.0/referLength;
    }
  
    //-----------------------------------------------------------------------------------------------

    private static void read_files(AtomList reference, AtomList protein1) {
    	AtomList al1,al2;
    	int ind1,ind2;
    	
    	len = 0;
    	al1 = reference.CAFilter();
    	al2 = protein1.CAFilter();
        for (int i=0 ; i<al1.size() ; i++) {
        	ind1 = al1.atomAt(i).residueNumber();
   	        for (int j=0 ; j<al2.size() ; j++) {
   	           ind2 = al2.atomAt(j).residueNumber();
   	           if (ind1==ind2) {
   	           	  len++;
   	           	  if (!al2.atomAt(j).residueName().equals(al1.atomAt(i).residueName())) {
   	           	     System.out.print(al2.atomAt(j) + "\n" + al1.atomAt(i)+ "\n");
   	           	     throw new RuntimeException("The residue names in both proteins mismatch:" +
   	           	     ind1 + " and " + ind2 + "\n\n");
   	           	  }
   	           }
   	        }
   	    }
   	    C1 = new double[3][len];
   	    C2 = new double[3][len];
   	    len = 0;
        for (int i=0 ; i<al1.size() ; i++) {
        	ind1 = al1.atomAt(i).residueNumber();
   	        for (int j=0 ; j<al2.size() ; j++) {
   	           ind2 = al2.atomAt(j).residueNumber();
   	           if (ind1==ind2) {
   	              C1[0][len] =  al1.atomAt(i).x();
   	              C1[1][len] =  al1.atomAt(i).y();
   	              C1[2][len] =  al1.atomAt(i).z();
   	              C2[0][len] =  al2.atomAt(j).x();
   	              C2[1][len] =  al2.atomAt(j).y();
   	              C2[2][len] =  al2.atomAt(j).z();
   	           	  len++; 
   	           	}
   	        }
   	    }
    }//read_files
    
    
}//class 
    
    
	



        	
   
