package programs;

import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.applications.homologyModelingNir.ProbA;
import meshi.applications.homologyModelingNir.TrivialHomologyModeling;
import meshi.applications.prediction.GDTcalculator;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.ROT1solvation.CBSolvationCreator;
import meshi.energy.ROT1solvation.CentroidSolvationCreator;
import meshi.energy.ROT1solvation.ROT1SolvationCreator;
import meshi.energy.ROT1solvation.ROT1SolvationEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.ResidueBuilder;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.AtomTypes;
import meshi.parameters.MeshiPotential;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;
import meshi.util.file.MeshiWriter;
import meshi.util.rotamericTools.RotamericTools;

public class AnalyzingSuboptimalsAllResidues extends MeshiProgram implements Residues, 
	AtomTypes , KeyWords, MeshiPotential { 
	
	
	static String commandsFileName;
	static String dirName;

		public static void main(String[] args)  {
			init(args); 
			CommandList commands = new CommandList(commandsFileName);
			ProbA probA = new ProbA(dirName + "/probA_output");
			String[] probA_input = File2StringArray.f2a(dirName + "/probA_input");
			int templateStart = (new Integer(probA_input[5])).intValue();
			Protein solution = new Protein(new AtomList(dirName+"/Native.pdb"),new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
			Protein template = new Protein(new AtomList(dirName+"/Template.pdb"),new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
			DunbrackLib lib = new DunbrackLib(commands, 1.0, 1);
			double[][] ppTemplate = RotamericTools.phipsi(template, new DistanceMatrix(template.atoms(), 4, 4.0, 2.0));
			
			String[] psiblastAlignment = File2StringArray.f2a(dirName+"/psiblast");
			Protein psiblastHomology = TrivialHomologyModeling.trivialHomology(template, templateStart, 
					psiblastAlignment[0], psiblastAlignment[1], 
					ppTemplate, lib);
			double psiblastGDT = GDTcalculator.gdt(solution.atoms(), psiblastHomology.atoms(),  1.0, 2.0, 4.0, 8.0);			
			System.out.println("PsiBlast GDT: " + psiblastGDT);
			double psiblastGDTHA = GDTcalculator.gdt(solution.atoms(), psiblastHomology.atoms(), 0.5,  1.0, 2.0, 4.0);			
			System.out.println("PsiBlast GDTHA: " + psiblastGDTHA);
			double psiblastGDT3 = GDTcalculator.gdt(solution.atoms(), psiblastHomology.atoms(), 3.0,  3.0, 3.0, 3.0);			
			System.out.println("PsiBlast GDT3: " + psiblastGDT3);
			try {
			psiblastHomology.atoms().print(new MeshiWriter(dirName+"/psiHomology.pdb"));
		} catch (Exception e) {			}
			
			String[] structuralAlignment = File2StringArray.f2a(dirName+"/structural");
			Protein structuralHomology = TrivialHomologyModeling.trivialHomology(template, templateStart, 
					structuralAlignment[0], structuralAlignment[1], 
					ppTemplate, lib);
			double structuralGDT = GDTcalculator.gdt(solution.atoms(), structuralHomology.atoms(),  1.0, 2.0, 4.0, 8.0);			
			System.out.println("structural GDT: " + structuralGDT);
			double structuralGDTHA = GDTcalculator.gdt(solution.atoms(), structuralHomology.atoms(), 0.5,  1.0, 2.0, 4.0);			
			System.out.println("structural GDTHA: " + structuralGDTHA);
			double structuralGDT3 = GDTcalculator.gdt(solution.atoms(), structuralHomology.atoms(),  3.0, 3.0, 3.0, 3.0);			
			System.out.println("structural GDT3: " + structuralGDT3);
			try {
				structuralHomology.atoms().print(new MeshiWriter(dirName+"/strHomology.pdb"));
			} catch (Exception e) {			}
			
			double[][] gdts = new double[probA.alignmentNum()+2][3];
			double[] E_ROT1_baysian_50 = new double[probA.alignmentNum()+2];
			double[] E_ROT1_baysian_60 = new double[probA.alignmentNum()+2];
			double[] E_ROT1_baysian_70 = new double[probA.alignmentNum()+2];
			double[] E_ROT1_baysian_80 = new double[probA.alignmentNum()+2];
			double[] E_CB_baysian_50 = new double[probA.alignmentNum()+2];
			double[] E_CB_baysian_60 = new double[probA.alignmentNum()+2];
			double[] E_CB_baysian_70 = new double[probA.alignmentNum()+2];
			double[] E_CB_baysian_80 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_baysian_50 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_baysian_60 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_baysian_70 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_baysian_80 = new double[probA.alignmentNum()+2];
			double[] E_ROT1_SA_50 = new double[probA.alignmentNum()+2];
			double[] E_ROT1_SA_60 = new double[probA.alignmentNum()+2];
			double[] E_ROT1_SA_70 = new double[probA.alignmentNum()+2];
			double[] E_ROT1_SA_80 = new double[probA.alignmentNum()+2];
			double[] E_CB_SA_50 = new double[probA.alignmentNum()+2];
			double[] E_CB_SA_60 = new double[probA.alignmentNum()+2];
			double[] E_CB_SA_70 = new double[probA.alignmentNum()+2];
			double[] E_CB_SA_80 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_SA_50 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_SA_60 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_SA_70 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_SA_80 = new double[probA.alignmentNum()+2];
			double[] E_ROT1_SAphob_50 = new double[probA.alignmentNum()+2];
			double[] E_ROT1_SAphob_60 = new double[probA.alignmentNum()+2];
			double[] E_ROT1_SAphob_70 = new double[probA.alignmentNum()+2];
			double[] E_ROT1_SAphob_80 = new double[probA.alignmentNum()+2];
			double[] E_CB_SAphob_50 = new double[probA.alignmentNum()+2];
			double[] E_CB_SAphob_60 = new double[probA.alignmentNum()+2];
			double[] E_CB_SAphob_70 = new double[probA.alignmentNum()+2];
			double[] E_CB_SAphob_80 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_SAphob_50 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_SAphob_60 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_SAphob_70 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_SAphob_80 = new double[probA.alignmentNum()+2];
			double[] E_ROT1_REP_50 = new double[probA.alignmentNum()+2];
			double[] E_ROT1_REP_60 = new double[probA.alignmentNum()+2];
			double[] E_ROT1_REP_70 = new double[probA.alignmentNum()+2];
			double[] E_ROT1_REP_80 = new double[probA.alignmentNum()+2];
			double[] E_CB_REP_50 = new double[probA.alignmentNum()+2];
			double[] E_CB_REP_60 = new double[probA.alignmentNum()+2];
			double[] E_CB_REP_70 = new double[probA.alignmentNum()+2];
			double[] E_CB_REP_80 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_REP_50 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_REP_60 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_REP_70 = new double[probA.alignmentNum()+2];
			double[] E_CENTROID_REP_80 = new double[probA.alignmentNum()+2];

			ROT1SolvationEnergy rot1_50 = null;
			ROT1SolvationEnergy rot1_60 = null;
			ROT1SolvationEnergy rot1_70 = null;
			ROT1SolvationEnergy rot1_80 = null;
			ROT1SolvationEnergy cb_50 = null;
			ROT1SolvationEnergy cb_60 = null;
			ROT1SolvationEnergy cb_70 = null;
			ROT1SolvationEnergy cb_80 = null;
			ROT1SolvationEnergy centroid_50 = null;
			ROT1SolvationEnergy centroid_60 = null;
			ROT1SolvationEnergy centroid_70 = null;
			ROT1SolvationEnergy centroid_80 = null;

			EnergyCreator[] energyCreators = {new ROT1SolvationCreator(1.0,"5.0",false),
					new ROT1SolvationCreator(1.0,"6.0",false),
					new ROT1SolvationCreator(1.0,"7.0",false),
					new ROT1SolvationCreator(1.0,"8.0",false),
					new CBSolvationCreator(1.0,"5.0",false),
					new CBSolvationCreator(1.0,"6.0",false),
					new CBSolvationCreator(1.0,"7.0",false),
					new CBSolvationCreator(1.0,"8.0",false),
					new CentroidSolvationCreator(1.0,"5.0",false),
					new CentroidSolvationCreator(1.0,"6.0",false),
					new CentroidSolvationCreator(1.0,"7.0",false),
					new CentroidSolvationCreator(1.0,"8.0",false)};
			
			for (int ali=0 ; ali<probA.alignmentNum()+2 ; ali++) {
				System.out.println("******************\n" + ali + "\n******************\n");
				Protein homology;
				if (ali==probA.alignmentNum())
					homology = psiblastHomology;
				else if (ali==(probA.alignmentNum()+1))
					homology = structuralHomology;
				else
					homology = TrivialHomologyModeling.trivialHomology(template, templateStart, 
						probA.suboptimalTops(ali), probA.suboptimalBottoms(ali), 
						ppTemplate, lib);
				gdts[ali][0] = GDTcalculator.gdt(solution.atoms(), homology.atoms(),  1.0, 2.0, 4.0, 8.0);
				gdts[ali][1] = GDTcalculator.gdt(solution.atoms(), homology.atoms(),  0.5 ,1.0, 2.0, 4.0);
				gdts[ali][2] = GDTcalculator.gdt(solution.atoms(), homology.atoms(),  3.0, 3.0, 3.0, 3.0);
				TotalEnergy energy = new TotalEnergy(homology, new DistanceMatrix(homology.atoms(),  8.0, 0.1, 4), energyCreators, commands);
				rot1_50 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[0]);
				rot1_60 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[1]);
				rot1_70 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[2]);
				rot1_80 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[3]);
				cb_50 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[4]);
				cb_60 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[5]);
				cb_70 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[6]);
				cb_80 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[7]);
				centroid_50 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[8]);
				centroid_60 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[9]);
				centroid_70 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[10]);
				centroid_80 = (ROT1SolvationEnergy) (energy.getEnergyTerms(new ROT1SolvationEnergy())[11]);

				E_ROT1_baysian_50[ali] = rot1_50.evaluate();
				E_ROT1_baysian_60[ali] = rot1_60.evaluate();
				E_ROT1_baysian_70[ali] = rot1_70.evaluate();
				E_ROT1_baysian_80[ali] = rot1_80.evaluate();
				E_CB_baysian_50[ali] = cb_50.evaluate();
				E_CB_baysian_60[ali] = cb_60.evaluate();
				E_CB_baysian_70[ali] = cb_70.evaluate();
				E_CB_baysian_80[ali] = cb_80.evaluate();

				E_ROT1_SAphob_50[ali] = rot1_50.getEnergyHydrophobicWeightedSA();
				E_ROT1_SAphob_60[ali] = rot1_60.getEnergyHydrophobicWeightedSA();
				E_ROT1_SAphob_70[ali] = rot1_70.getEnergyHydrophobicWeightedSA();
				E_ROT1_SAphob_80[ali] = rot1_80.getEnergyHydrophobicWeightedSA();
				E_CB_SAphob_50[ali] = cb_50.getEnergyHydrophobicWeightedSA();
				E_CB_SAphob_60[ali] = cb_60.getEnergyHydrophobicWeightedSA();
				E_CB_SAphob_70[ali] = cb_70.getEnergyHydrophobicWeightedSA();
				E_CB_SAphob_80[ali] = cb_80.getEnergyHydrophobicWeightedSA();

				E_ROT1_SA_50[ali] = E_ROT1_SAphob_50[ali] + rot1_50.getEnergyPolarWeightedSA();
				E_ROT1_SA_60[ali] = E_ROT1_SAphob_60[ali] + rot1_60.getEnergyPolarWeightedSA();
				E_ROT1_SA_70[ali] = E_ROT1_SAphob_70[ali] + rot1_70.getEnergyPolarWeightedSA();
				E_ROT1_SA_80[ali] = E_ROT1_SAphob_80[ali] + rot1_80.getEnergyPolarWeightedSA();
				E_CB_SA_50[ali] = E_CB_SAphob_50[ali] + cb_50.getEnergyPolarWeightedSA();
				E_CB_SA_60[ali] = E_CB_SAphob_60[ali] + cb_60.getEnergyPolarWeightedSA();
				E_CB_SA_70[ali] = E_CB_SAphob_70[ali] + cb_70.getEnergyPolarWeightedSA();
				E_CB_SA_80[ali] = E_CB_SAphob_80[ali] + cb_80.getEnergyPolarWeightedSA();
				
				E_ROT1_REP_50[ali] = rot1_50.getEnergyRepresentativePolar();
				E_ROT1_REP_60[ali] = rot1_60.getEnergyRepresentativePolar();
				E_ROT1_REP_70[ali] = rot1_70.getEnergyRepresentativePolar();
				E_ROT1_REP_80[ali] = rot1_80.getEnergyRepresentativePolar();
				E_CB_REP_50[ali] = cb_50.getEnergyRepresentativePolar();
				E_CB_REP_60[ali] = cb_60.getEnergyRepresentativePolar();
				E_CB_REP_70[ali] = cb_70.getEnergyRepresentativePolar();
				E_CB_REP_80[ali] = cb_80.getEnergyRepresentativePolar();
				
				
				for (int res=0 ; res<homology.residues().size(); res++) {
					if ((homology.residues().residueAt(res).type<20) &&
							(homology.residues().residueAt(res).type>-1))
						ResidueBuilder.buildCentroid(homology.residues().residueAt(res));
				}
				energy.evaluate();
				E_CENTROID_baysian_50[ali] = centroid_50.evaluate();
				E_CENTROID_baysian_60[ali] = centroid_60.evaluate();
				E_CENTROID_baysian_70[ali] = centroid_70.evaluate();
				E_CENTROID_baysian_80[ali] = centroid_80.evaluate();
				E_CENTROID_SAphob_50[ali] = centroid_50.getEnergyHydrophobicWeightedSA();
				E_CENTROID_SAphob_60[ali] = centroid_60.getEnergyHydrophobicWeightedSA();
				E_CENTROID_SAphob_70[ali] = centroid_70.getEnergyHydrophobicWeightedSA();
				E_CENTROID_SAphob_80[ali] = centroid_80.getEnergyHydrophobicWeightedSA();
				E_CENTROID_SA_50[ali] = E_CENTROID_SAphob_50[ali] + centroid_50.getEnergyPolarWeightedSA();
				E_CENTROID_SA_60[ali] = E_CENTROID_SAphob_60[ali] + centroid_60.getEnergyPolarWeightedSA();
				E_CENTROID_SA_70[ali] = E_CENTROID_SAphob_70[ali] + centroid_70.getEnergyPolarWeightedSA();
				E_CENTROID_SA_80[ali] = E_CENTROID_SAphob_80[ali] + centroid_80.getEnergyPolarWeightedSA();
				E_CENTROID_REP_50[ali] = centroid_50.getEnergyRepresentativePolar();
				E_CENTROID_REP_60[ali] = centroid_60.getEnergyRepresentativePolar();
				E_CENTROID_REP_70[ali] = centroid_70.getEnergyRepresentativePolar();
				E_CENTROID_REP_80[ali] = centroid_80.getEnergyRepresentativePolar();
				
			}

			
			
			for (int ali=0 ; ali<probA.alignmentNum()+2 ; ali++) {
				if (ali<probA.alignmentNum())
					System.out.println("99990000 " + ali + " " + gdts[ali][0] + " " + gdts[ali][1] + " " + gdts[ali][2] + " " + 
						probA.scores(ali));
				else
					System.out.println("99990000 " + ali + " " + gdts[ali][0] + " " + gdts[ali][1] + " " + gdts[ali][2] + " 0");					
				System.out.println("99990500 " + ali + " " + E_ROT1_baysian_50[ali] + " " + E_ROT1_SA_50[ali] + " " + E_ROT1_SAphob_50[ali] +
						" " + E_ROT1_REP_50[ali]);
				System.out.println("99990501 " + ali + " " + E_CB_baysian_50[ali] + " " + E_CB_SA_50[ali] + " " + E_CB_SAphob_50[ali] +
						" " + E_CB_REP_50[ali]);
				System.out.println("99990502 " + ali + " " + E_CENTROID_baysian_50[ali] + " " + E_CENTROID_SA_50[ali] + " " + E_CENTROID_SAphob_50[ali] +
						" " + E_CENTROID_REP_50[ali]);
				System.out.println("99990600 " + ali + " " + E_ROT1_baysian_60[ali] + " " + E_ROT1_SA_60[ali] + " " + E_ROT1_SAphob_60[ali] +
						" " + E_ROT1_REP_60[ali]);
				System.out.println("99990601 " + ali + " " + E_CB_baysian_60[ali] + " " + E_CB_SA_60[ali] + " " + E_CB_SAphob_60[ali] +
						" " + E_CB_REP_60[ali]);
				System.out.println("99990602 " + ali + " " + E_CENTROID_baysian_60[ali] + " " + E_CENTROID_SA_60[ali] + " " + E_CENTROID_SAphob_60[ali] +
						" " + E_CENTROID_REP_60[ali]);
				System.out.println("99990700 " + ali + " " + E_ROT1_baysian_70[ali] + " " + E_ROT1_SA_70[ali] + " " + E_ROT1_SAphob_70[ali] +
						" " + E_ROT1_REP_70[ali]);
				System.out.println("99990701 " + ali + " " + E_CB_baysian_70[ali] + " " + E_CB_SA_70[ali] + " " + E_CB_SAphob_70[ali] +
						" " + E_CB_REP_70[ali]);
				System.out.println("99990702 " + ali + " " + E_CENTROID_baysian_70[ali] + " " + E_CENTROID_SA_70[ali] + " " + E_CENTROID_SAphob_70[ali] +
						" " + E_CENTROID_REP_70[ali]);
				System.out.println("99990800 " + ali + " " + E_ROT1_baysian_80[ali] + " " + E_ROT1_SA_80[ali] + " " + E_ROT1_SAphob_80[ali] +
						" " + E_ROT1_REP_80[ali]);
				System.out.println("99990801 " + ali + " " + E_CB_baysian_80[ali] + " " + E_CB_SA_80[ali] + " " + E_CB_SAphob_80[ali] +
						" " + E_CB_REP_80[ali]);
				System.out.println("99990802 " + ali + " " + E_CENTROID_baysian_80[ali] + " " + E_CENTROID_SA_80[ali] + " " + E_CENTROID_SAphob_80[ali] +
						" " + E_CENTROID_REP_80[ali]);
			}
			
			try{
				BufferedWriter bw = new BufferedWriter(new FileWriter(dirName + "/energies.txt"));
				for (int ali=0 ; ali<probA.alignmentNum()+2 ; ali++) {
					if (ali<probA.alignmentNum())
						bw.write("99990000 " + ali + " " + gdts[ali][0] + " " + gdts[ali][1] + " " + gdts[ali][2] + " " + 
							probA.scores(ali) + "\n");
					else
						bw.write("99990000 " + ali + " " + gdts[ali][0] + " " + gdts[ali][1] + " " + gdts[ali][2] + " 0\n");					
					bw.write("99990500 " + ali + " " + E_ROT1_baysian_50[ali] + " " + E_ROT1_SA_50[ali] + " " + E_ROT1_SAphob_50[ali] +
							" " + E_ROT1_REP_50[ali] + "\n");
					bw.write("99990501 " + ali + " " + E_CB_baysian_50[ali] + " " + E_CB_SA_50[ali] + " " + E_CB_SAphob_50[ali] +
							" " + E_CB_REP_50[ali] + "\n");
					bw.write("99990502 " + ali + " " + E_CENTROID_baysian_50[ali] + " " + E_CENTROID_SA_50[ali] + " " + E_CENTROID_SAphob_50[ali] +
							" " + E_CENTROID_REP_50[ali] + "\n");
					bw.write("99990600 " + ali + " " + E_ROT1_baysian_60[ali] + " " + E_ROT1_SA_60[ali] + " " + E_ROT1_SAphob_60[ali] +
							" " + E_ROT1_REP_60[ali] + "\n");
					bw.write("99990601 " + ali + " " + E_CB_baysian_60[ali] + " " + E_CB_SA_60[ali] + " " + E_CB_SAphob_60[ali] +
							" " + E_CB_REP_60[ali] + "\n");
					bw.write("99990602 " + ali + " " + E_CENTROID_baysian_60[ali] + " " + E_CENTROID_SA_60[ali] + " " + E_CENTROID_SAphob_60[ali] +
							" " + E_CENTROID_REP_60[ali] + "\n");
					bw.write("99990700 " + ali + " " + E_ROT1_baysian_70[ali] + " " + E_ROT1_SA_70[ali] + " " + E_ROT1_SAphob_70[ali] +
							" " + E_ROT1_REP_70[ali] + "\n");
					bw.write("99990701 " + ali + " " + E_CB_baysian_70[ali] + " " + E_CB_SA_70[ali] + " " + E_CB_SAphob_70[ali] +
							" " + E_CB_REP_70[ali] + "\n");
					bw.write("99990702 " + ali + " " + E_CENTROID_baysian_70[ali] + " " + E_CENTROID_SA_70[ali] + " " + E_CENTROID_SAphob_70[ali] +
							" " + E_CENTROID_REP_70[ali] + "\n");
					bw.write("99990800 " + ali + " " + E_ROT1_baysian_80[ali] + " " + E_ROT1_SA_80[ali] + " " + E_ROT1_SAphob_80[ali] +
							" " + E_ROT1_REP_80[ali] + "\n");
					bw.write("99990801 " + ali + " " + E_CB_baysian_80[ali] + " " + E_CB_SA_80[ali] + " " + E_CB_SAphob_80[ali] +
							" " + E_CB_REP_80[ali] + "\n");
					bw.write("99990802 " + ali + " " + E_CENTROID_baysian_80[ali] + " " + E_CENTROID_SA_80[ali] + " " + E_CENTROID_SAphob_80[ali] +
							" " + E_CENTROID_REP_80[ali] + "\n");
				}
				bw.close();
			}
			catch(Exception e) {
				throw new RuntimeException(e.getMessage());
			}
	  
			//			try {
//				homology.atoms().print(new MeshiWriter(dirName+"/homology.pdb"));
//			} catch (IOException e) {			}
			
								
		}
			
			


		/** ================================= init =========================================
		 *
		 *A static function for parsing of the command line arguments and assigning the 
		 *variables commandsFileName, modelFileName and randomNumberSeed with the right inputs. Note that this
		 *static method is using parsing functions such as getOrderedArguments that are defined in MeshiProgram
		 *that MinimizeProtein inherits.
		 **/

		protected static void init(String[] args) {

			int zvl = ALA; // force the reading of "meshi.parameters.Residues"
			zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"


			String errorMessage = ("\n                  ******************\n"+
					"Usage java -Xmx300m AnalyzingSuboptimal <commands file name> \n"+
			"                    ******************\n");

			if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));

			commandsFileName = getOrderedArgument(args);
			if (commandsFileName == null) throw new RuntimeException(errorMessage);
			System.out.println("# commandsFileName = "+commandsFileName);
			
			dirName = getOrderedArgument(args);
			if (dirName == null) throw new RuntimeException(errorMessage);
			System.out.println("# Working on directory: "+dirName);			

			initRandom(0);

		}
	
}
