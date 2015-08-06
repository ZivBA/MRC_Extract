package meshi.util.TRIC.interRing;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import meshi.molecularElements.AtomList;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.TRIC.PutUnitInAnyTopPosition;

public class PrepareInterRingPair extends MeshiProgram implements Residues { 
   
	private static String pathToInterRingFolder = null;

	private static String unitLetterUP = null;
 
	private static String unitLetterDOWN = null;

	private static String unitLetterHomology = null;

    private static String outputFile = null;  
 
    public static void main(String[] args) {
	init(args); 

	char convertHomologyToThermosome = 0;
	if (unitLetterHomology.equals("A"))
		convertHomologyToThermosome = 'I';
	else if (unitLetterHomology.equals("B"))
		convertHomologyToThermosome = 'J';
	else 
		throw new RuntimeException("Only the letters A and B are allowed as thermosome homology units.");

	String writeDir = pathToInterRingFolder+"/inter_ring/"+unitLetterUP.charAt(0)+
	unitLetterDOWN.charAt(0);

	// Doing the UP unit
	AtomList statUP = new AtomList(pathToInterRingFolder+"/inter_ring/aux_files/up_"+unitLetterHomology.charAt(0)+".pdb");
	AtomList statDOWN = new AtomList(pathToInterRingFolder+"/inter_ring/aux_files/down_"+unitLetterHomology.charAt(0)+".pdb");
	AtomList moveUP = new AtomList(pathToInterRingFolder+"/inter_ring/aux_files/unit_"+unitLetterUP.charAt(0)+".pdb");
	AtomList moveDOWN = new AtomList(pathToInterRingFolder+"/inter_ring/aux_files/unit_"+unitLetterDOWN.charAt(0)+".pdb");
	PutUnitInAnyTopPosition.alignByEquatorialDomains(statUP, convertHomologyToThermosome,
			moveUP, unitLetterUP.charAt(0));
	PutUnitInAnyTopPosition.alignByEquatorialDomains(statDOWN, convertHomologyToThermosome,
			moveDOWN, unitLetterDOWN.charAt(0));
	for (int atomC=0 ; atomC<moveUP.size() ; atomC++) {
		moveUP.atomAt(atomC).setChain("A");
	}
	for (int atomC=0 ; atomC<moveDOWN.size() ; atomC++) {
		moveDOWN.atomAt(atomC).setChain("B");
	}
		
	try{
		// Creating the directory and writing to disk
		boolean success = (new File(writeDir)).mkdir();
		if (success) {
			System.out.println("Directory: " + writeDir + " created");
		}    
		BufferedWriter bw = new BufferedWriter(new FileWriter(writeDir + "/" + outputFile + ".pdb"));
		for (int atomC=0 ; atomC<moveUP.size() ; atomC++) {
			bw.write(moveUP.atomAt(atomC) + "\n");
		}
		bw.write("TER\n");
		for (int atomC=0 ; atomC<moveDOWN.size() ; atomC++) {
			bw.write(moveDOWN.atomAt(atomC) + "\n");
		}
		bw.write("TER\n");
		bw.write("END\n");
		bw.close();
	}
	catch(Exception e) {
	    throw new RuntimeException(e.getMessage());
	}    	

	// Creating a separated complex
	for (int atomC=0 ; atomC<moveUP.size() ; atomC++) {
		moveUP.atomAt(atomC).setZ(moveUP.atomAt(atomC).z()-20.0);
	}

	try{
		// Writing the separated complex to disk
		BufferedWriter bw = new BufferedWriter(new FileWriter(writeDir + "/" + outputFile + ".far_apart.pdb"));
		for (int atomC=0 ; atomC<moveUP.size() ; atomC++) {
			bw.write(moveUP.atomAt(atomC) + "\n");
		}
		bw.write("TER\n");
		for (int atomC=0 ; atomC<moveDOWN.size() ; atomC++) {
			bw.write(moveDOWN.atomAt(atomC) + "\n");
		}
		bw.write("TER\n");
		bw.write("END\n");
		bw.close();
	}
	catch(Exception e) {
	    throw new RuntimeException(e.getMessage());
	}    		
}
    

    /** ================================= init =========================================
     *
     *A static function for parsing of the command line arguments and assigning the 
     *variables commandsFileName, modelFileName and randomNumberSeed with the right inputs. Note that this
     *static method is using parsing functions such as getOrderedArguments that are defined in MeshiProgram
     *that MinimizeProtein inherits.
     **/
     
    protected static void init(String[] args) {
 
	/**** NOTE *** the next two lines. Because of a BUG in the Java VM, the 
	 * interfaces "Residues" and "AtomTypes" are not loaded automatically when MinimizeProtein initialize. 
	 * For this purpose these two lines are crucial wherever these two interfaces are implemented. The user might 
	 * rightfully feel that these two lines are "black magic" programming, but happily to our knowledge this is 
	 * the only bizarre phenomenon we are aware of in meshi.
	 **/

	String errorMessage = ("\n                  ******************\n"+
			       "Usage java  PrepareInterRingPair <Path to the inter_ring folder (without the inter_ring in the end)> " +
			       "<Unit letter for UP: AZBGQDEHIJ> <Unit letter for DOWN: AZBGQDEHIJ> <which thermosome homologue: A,B> " +
			       "<output file prefix for the complex> \n"+
			       "                    ******************\n");
			      
	pathToInterRingFolder = getOrderedArgument(args);
	if (pathToInterRingFolder == null) throw new RuntimeException(errorMessage);
	pathToInterRingFolder = pathToInterRingFolder.trim();
	System.out.println("# The path to the inter_ring folder: " + pathToInterRingFolder);

	unitLetterUP = getOrderedArgument(args);
	if (unitLetterUP == null) throw new RuntimeException(errorMessage);
	unitLetterUP = unitLetterUP.trim();
	System.out.println("# The UP unit is: " + unitLetterUP);

	unitLetterDOWN = getOrderedArgument(args);
	if (unitLetterDOWN == null) throw new RuntimeException(errorMessage);
	unitLetterDOWN = unitLetterDOWN.trim();
	System.out.println("# The DOWN unit is: " + unitLetterDOWN);

	unitLetterHomology = getOrderedArgument(args);
	if (unitLetterHomology == null) throw new RuntimeException(errorMessage);
	unitLetterHomology = unitLetterHomology.trim();
	System.out.println("# The homology is based on thermosome unit: " + unitLetterHomology);

	outputFile = getOrderedArgument(args);
	if (outputFile == null) throw new RuntimeException(errorMessage);
	System.out.println("# Output file name is "+outputFile);
	
	initRandom(999);
    }
}
