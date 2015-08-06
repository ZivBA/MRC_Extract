package meshi.util.TRIC.interRing;

import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;

public class RandomRotTrans extends MeshiProgram implements Residues { 
   
	private static String initialStruct = null;

	private static double translateRange = 0.0;
 
	private static double translateJumps = 0.0;

	private static double rotateRange = 0.0;
	 
	private static double rotateJumps = 0.0;

	private static String outputPrefix = null;
 
    public static void main(String[] args) {
	init(args); 

	int structureCounter = 0;
	AtomList complexOrig = new AtomList(initialStruct);
	AtomList unitAorig = complexOrig.filter(new AtomList.ChainA());
	AtomList unitBorig = complexOrig.filter(new AtomList.ChainB());
	AtomList complexMove = new AtomList(initialStruct);
	AtomList unitBmove = complexMove.filter(new AtomList.ChainB());
	// Calculating CM for unit B
	double cmx, cmy, cmz; // center of mass x, y and z
	cmx = cmy = cmz = 0.0;
	for (int c=0 ; c<unitBorig.size() ; c++) {
		cmx += unitBorig.atomAt(c).x();
		cmy += unitBorig.atomAt(c).y();
		cmz += unitBorig.atomAt(c).z();
	}
	cmx /= unitBorig.size();
	cmy /= unitBorig.size();
	cmz /= unitBorig.size();
	
	printToFile(structureCounter,unitAorig,unitBmove);
	
	for (double moveX = -translateRange ; moveX<(translateRange+1e-6) ; moveX+=translateJumps) {
		for (double moveY = -translateRange ; moveY<(translateRange+1e-6) ; moveY+=translateJumps) {
			for (double moveZ = -0.25 ; moveZ<(0.25+1e-6) ; moveZ+=0.25) {
				for (double rotX = -rotateRange ; rotX<(rotateRange+1e-6) ; rotX+=rotateJumps) {
					for (double rotY = -rotateRange ; rotY<(rotateRange+1e-6) ; rotY+=rotateJumps) {
						for (double rotZ = -rotateRange ; rotZ<(rotateRange+1e-6) ; rotZ+=rotateJumps) {
							if (Math.sqrt(rotX*rotX+rotY*rotY+rotZ*rotZ) < (1.42*rotateRange)) {
								if (Math.sqrt(moveX*moveX+moveY*moveY+moveZ*moveZ) < (1.32*translateRange)) {
									structureCounter++;
									System.out.println(moveX+" "+moveY+" "+moveZ+" "+rotX+" "+rotY+" "+rotZ+" "+structureCounter);
									transAndRot(unitBmove, cmx, cmy, cmz, moveX, moveY, moveZ, rotX, rotY, rotZ);
									printToFile(structureCounter,unitAorig,unitBmove);
									for (int c=0 ; c<unitBmove.size() ; c++) {
										Atom atom = unitBorig.findAtomInList(unitBmove.atomAt(c).name(), unitBmove.atomAt(c).residueNumber());
										unitBmove.atomAt(c).setXYZ(atom.x(), atom.y(), atom.z());								
									}
								}
							}
						}
					}					
				}				
			}			
		}
	}
	
}
    
    
    
    static void transAndRot(AtomList unitB, double cmx, double cmy, double cmz, 
    		double moveX, double moveY, double moveZ, double rotX, double rotY, double rotZ) {
    	for (int c=0 ; c<unitB.size() ; c++) {
    		unitB.atomAt(c).setX(unitB.atomAt(c).x() - cmx);
    		unitB.atomAt(c).setY(unitB.atomAt(c).y() - cmy);
    		unitB.atomAt(c).setZ(unitB.atomAt(c).z() - cmz);
    		double tmpx = 1.0*unitB.atomAt(c).x();
    		double tmpy = Math.cos(rotX)*unitB.atomAt(c).y() + Math.sin(rotX)*unitB.atomAt(c).z();
    		double tmpz = Math.sin(-rotX)*unitB.atomAt(c).y() + Math.cos(rotX)*unitB.atomAt(c).z();
    		unitB.atomAt(c).setX(tmpx);
    		unitB.atomAt(c).setY(tmpy);
    		unitB.atomAt(c).setZ(tmpz);
    		tmpx = Math.cos(rotY)*unitB.atomAt(c).x() + Math.sin(rotY)*unitB.atomAt(c).z();
    		tmpy = 1.0*unitB.atomAt(c).y();
    		tmpz = Math.sin(-rotY)*unitB.atomAt(c).x() + Math.cos(rotY)*unitB.atomAt(c).z();
    		unitB.atomAt(c).setX(tmpx);
    		unitB.atomAt(c).setY(tmpy);
    		unitB.atomAt(c).setZ(tmpz);
    		tmpx = Math.cos(rotZ)*unitB.atomAt(c).x() + Math.sin(rotZ)*unitB.atomAt(c).y();
    		tmpy = Math.sin(-rotZ)*unitB.atomAt(c).x() + Math.cos(rotZ)*unitB.atomAt(c).y();
    		tmpz = 1.0*unitB.atomAt(c).z();
    		unitB.atomAt(c).setX(tmpx);
    		unitB.atomAt(c).setY(tmpy);
    		unitB.atomAt(c).setZ(tmpz);
    		unitB.atomAt(c).setX(unitB.atomAt(c).x() + cmx + moveX);
    		unitB.atomAt(c).setY(unitB.atomAt(c).y() + cmy + moveY);
    		unitB.atomAt(c).setZ(unitB.atomAt(c).z() + cmz + moveZ);    		
    	}    	
    }
    
	static void printToFile(int structureCounter, AtomList unitA, AtomList unitB) {
		try{
			// Writing the separated complex to disk
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputPrefix + "." + structureCounter + ".pdb"));
			for (int atomC=0 ; atomC<unitA.size() ; atomC++) {
				bw.write(unitA.atomAt(atomC) + "\n");
			}
			bw.write("TER\n");
			for (int atomC=0 ; atomC<unitB.size() ; atomC++) {
				bw.write(unitB.atomAt(atomC) + "\n");
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
			       "Usage java  RandomRotTrans <Initial PDB> " +
			       "<Translate range> <Translate jumps> <Rotate range> <Rotate jumps> " +
			       "<output prefix (will add .XXX.pdb to end)> \n"+
			       "                    ******************\n");
			      
	initialStruct = getOrderedArgument(args);
	if (initialStruct == null) throw new RuntimeException(errorMessage);
	initialStruct = initialStruct.trim();
	System.out.println("# The initial structure is: " + initialStruct);

	String tmp = getOrderedArgument(args);
	if (tmp == null) throw new RuntimeException(errorMessage);
	translateRange = Double.parseDouble(tmp);
	System.out.println("# Translate range: " + translateRange);
	
	tmp = getOrderedArgument(args);
	if (tmp == null) throw new RuntimeException(errorMessage);
	translateJumps = Double.parseDouble(tmp);
	System.out.println("# Translate jump: " + translateJumps);

	tmp = getOrderedArgument(args);
	if (tmp == null) throw new RuntimeException(errorMessage);
	rotateRange = Double.parseDouble(tmp);
	System.out.println("# Rotate range: " + rotateRange);
	rotateRange *= Math.PI/180.0;	
	
	tmp = getOrderedArgument(args);
	if (tmp == null) throw new RuntimeException(errorMessage);
	rotateJumps = Double.parseDouble(tmp);
	System.out.println("# Rotate jump: " + rotateJumps);
	rotateJumps *= Math.PI/180.0;

	outputPrefix = getOrderedArgument(args);
	if (outputPrefix == null) throw new RuntimeException(errorMessage);
	System.out.println("# Output prefix is: " + outputPrefix);
	
	initRandom(999);
    }
}
