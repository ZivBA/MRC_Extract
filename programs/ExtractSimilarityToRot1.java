package programs;

import java.text.DecimalFormat;

import meshi.geometry.DistanceMatrix;
import meshi.geometry.ResidueBuilder;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;
import meshi.util.rotamericTools.RotamericTools;

public class ExtractSimilarityToRot1 extends MeshiProgram implements Residues, AtomTypes{ 

	// Parameters: 
	// ------------
	// The structure database:
	static String listOfStructures = "C:/Users/Nir/Loop_Building_Project/listPISCES.txt";

	public static void main(String[] args){
		init(args); 
		
		// Dunbrack Lib
		CommandList commands = new CommandList("commands");
		DunbrackLib lib = new DunbrackLib(commands, 1.0, 2);
		
		// Data arrays
		double[] numberOfRes = new double[20];
		double[] rot1Prob = new double[20];
		double[] under05 = new double[20];
		double[] under10 = new double[20];
		double[] under15 = new double[20];
		double[] under20 = new double[20];
		double[] under25 = new double[20];
		double[] under30 = new double[20];
		double[] under35 = new double[20];
		int totalResidues = 0;
		

		// Going over the models
		String[] models = File2StringArray.f2a( listOfStructures);
		for (int i=0 ; i<850 ; i++) { 		
			System.out.println("Reading: " + i + " " + models[i]);

			// Creating the model
			Protein model = new Protein(new AtomList("Pisces/" +models[i]), new ResidueExtendedAtoms(ADD_HYDROGENS_AND_FREEZE));
			model.defrost();
			double[][] pp = RotamericTools.phipsi(model, new DistanceMatrix(model.atoms(), 5.5 , 1.0, 4));
			try {
			for (int res=0 ; res<model.residues().size(); res++) {
				Residue residue = model.residues().residueAt(res);
				if ((residue.type<20) && (residue.type>0) && (residue.type!=5) &&
					RotamericTools.isSidechainComplete(residue)) {
						double[] repCoorsOrigin = RotamericTools.findRepCoors(residue);
						ResidueBuilder.build(model.residues().residueAt(res),
								model.residues().residueAt(res).type, 
								lib.getRotamer(model.residues().residueAt(res).type,
										pp[model.residues().residueAt(res).number][0], 
										pp[model.residues().residueAt(res).number][1], 0));
						double[] repCoorsROT1 = RotamericTools.findRepCoors(residue);
						numberOfRes[residue.type]++;
						rot1Prob[residue.type] += lib.getRotamerProb(model.residues().residueAt(res).type,
								pp[model.residues().residueAt(res).number][0], 
								pp[model.residues().residueAt(res).number][1], 0);
						double dis = Math.sqrt((repCoorsOrigin[0]-repCoorsROT1[0])*(repCoorsOrigin[0]-repCoorsROT1[0]) + 
								(repCoorsOrigin[1]-repCoorsROT1[1])*(repCoorsOrigin[1]-repCoorsROT1[1]) + 
								(repCoorsOrigin[2]-repCoorsROT1[2])*(repCoorsOrigin[2]-repCoorsROT1[2]));
						if (dis<=0.5)
							under05[residue.type]++;
						if (dis<=1.0)
							under10[residue.type]++;
						if (dis<=1.5)
							under15[residue.type]++;
						if (dis<=2.0)
							under20[residue.type]++;
						if (dis<=2.5)
							under25[residue.type]++;
						if (dis<=3.0)
							under30[residue.type]++;
						if (dis<=3.5)
							under35[residue.type]++;
						totalResidues++;
				}
			}
			}
			catch (Exception e) {
				System.out.println("Error in the processing of this protein. Continueing....");
			}
		}
		
		// Printing the results
		DecimalFormat fmt = new DecimalFormat("0.#");
		for (int c=0 ; c<20 ; c++) {
			if (numberOfRes[c]>0) {
				System.out.println(Residue.nameThreeLetters(c) + "  " + fmt.format(100.0*numberOfRes[c]/totalResidues) + " " + 
						fmt.format(100.0*rot1Prob[c]/numberOfRes[c]) + " " + fmt.format(100.0*under05[c]/numberOfRes[c]) + " " +
						fmt.format(100.0*under10[c]/numberOfRes[c]) + " " + fmt.format(100.0*under15[c]/numberOfRes[c]) + " " + 
						fmt.format(100.0*under20[c]/numberOfRes[c]) + " " + fmt.format(100.0*under25[c]/numberOfRes[c]) + " " + 
						fmt.format(100.0*under30[c]/numberOfRes[c]) + " " + fmt.format(100.0*under35[c]/numberOfRes[c]));
			}
		}
		System.out.println(totalResidues);

	}




	protected static void init(String[] args) {

		/**** NOTE *** the next two lines. Because of a BUG in the Java VM, the 
		 * interfaces "Residues" and "AtomTypes" are not loaded automatically when MinimizeProtein initialize. 
		 * For this purpose these two lines are crucial wherever these two interfaces are implemented. The user might 
		 * rightfully feel that these two lines are "black magic" programming, but happily to our knowledge this is 
		 * the only bizarre phenomenon we are aware of in meshi.
		 **/
		int zvl = ALA; // force the reading of "meshi.parameters.Residues"
		zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"

		String errorMessage = ("\n                  ******************\n"+
				"Usage java -Xmx1000m ExtractSimilarityToRot1\n"+
		"                    ******************\n");

		initRandom(333);
	}
}
