package meshi.applications.loopBuilding.applications;

import java.util.Vector;

import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.file.MeshiWriter;

public class ChopModelInPlaces extends MeshiProgram implements Residues, AtomTypes {


	public static void main(String[] args) {
		init(args); 
		
		String outputFileName = null;
		Protein output=null;
		Vector<Protein> loops = new Vector<Protein>();
		Vector<Integer> froms = new Vector<Integer>();
		Vector<Integer> tos = new Vector<Integer>();
		
		// Reading the command line
		try {
			outputFileName = args[0].trim();
			System.out.println("# Output to: "+outputFileName);
			int argCounter = 1;
			while (argCounter<args.length) {
				loops.add(new ExtendedAtomsProtein(args[argCounter].trim(),DO_NOT_ADD_ATOMS));
				System.out.println("# Loop " + ((argCounter-1)/3+1) + " model is: "+args[argCounter].trim());				
				froms.add(new Integer(args[argCounter+1].trim()));
				System.out.println("# Loop " + ((argCounter-1)/3+1) + " is from residue: "+args[argCounter+1].trim());
				tos.add(new Integer(args[argCounter+2].trim()));
				System.out.println("# Loop " + ((argCounter-1)/3+1) + " is till residue: "+args[argCounter+2].trim());
				argCounter += 3;
			}
		}
		catch (Exception e) {
			String errorMessage = ("\n                  ******************\n"+
					"Usage: java -Xmx300m ChopModelInPlaces <output filename> <loop #1 pdb> <loop #1 from res.> <loop #1 to res.> <loop #2 pdb> <loop #2 from res.> <loop #2 to res.>... \n"+
			"                    ******************\n");
			throw new RuntimeException(errorMessage);
		}
		
		// First building all the atom instances
		AtomList allAtoms = new AtomList();
		for (int loop=0 ; loop<loops.size() ; loop++){
			for (int res=froms.get(loop).intValue() ; res<=tos.get(loop).intValue() ; res++) {
				if ((loops.get(loop).residue(res)!=null) && (loops.get(loop).residue(res).ca()!=null)) {
					for (int atomC=0 ; atomC<loops.get(loop).residue(res).atoms().size() ; atomC++) {
						Atom atom = loops.get(loop).residue(res).atoms().atomAt(atomC);
						if (allAtoms.findAtomInList(atom.name(), res)==null) {
							allAtoms.add(new Atom(atom));
							System.out.println("Added: " + atom);
						}
						else
							allAtoms.findAtomInList(atom.name(), res).setXYZ(atom.x(), atom.y(), atom.z());
					}
				}
				else {
					System.out.println("WARNING: could not find residue " + res + " in loop number " + (loop+1));
				}
			}
		}
		// I cannot believe I am doing it, but we need to sort this list!
		
		output = new Protein(allAtoms.sortByResidueNumber(), new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		for (int cc=0 ; cc<output.atoms().size() ; cc++)
			output.atoms().atomAt(cc).setChain("A");

		// Writing to file
		try {
			output.atoms().print(new MeshiWriter(outputFileName));
		}
		catch (Exception e) {
			System.out.print("\nThere was a problem writing the output:\n" + e + "\n\nContinueing...\n\n");
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
		int zvl = ALA; // force the reading of "meshi.parameters.Residues"
		zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"



		initRandom(999);
	}	

} // Of AddMissingResidues
