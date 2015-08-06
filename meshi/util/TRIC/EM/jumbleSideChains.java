package meshi.util.TRIC.EM;

import java.io.IOException;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.compositeTorsions.compositePropensity2DwithPP.CompositePropensity2DCreator;
import meshi.energy.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.inflate.InflateCreator;
import meshi.energy.linearRG.LinearRgCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.softExcludedVol.SoftExcludedVol;
import meshi.energy.softExcludedVol.SoftExcludedVolCreator;
import meshi.energy.solvate.SolvateCreatorHBforMinimization;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.ResidueBuilder;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LBFGS;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.Minimizer;
import meshi.optimizers.MinimizerException;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.file.MeshiWriter;
import meshi.util.rotamericTools.RotamericTools;

public class jumbleSideChains extends MeshiProgram implements Residues {
	
	public static void main(String[] args) {
		initRandom();
		String[] names = {"Atemp_HM_A",
				"Atemp_HM_B",
				"Atemp_HM_G",
				"Atemp_HM_D",
				"Atemp_HM_E",
				"Atemp_HM_H",
				"Atemp_HM_Q",
				"Atemp_HM_Z"};
		for (int doProt=0; doProt<names.length ; doProt++) {
			Protein prot = new Protein(new AtomList(names[doProt]+".pdb"),new ResidueExtendedAtoms(ADD_ATOMS));
			RotamericTools.jumble(prot);
			EnergyCreator[] energyCreators = {  
					new BondCreator(),
					new AngleCreator(),
					new PlaneCreator(),
					new OutOfPlaneCreator(),
					new SoftExcludedVolCreator(100.0 , 4 , 1.0),
			};
			CommandList commands = new CommandList("C:\\Users\\Nir\\Loop_Building_Project\\commands");
			prot.atoms().backbone().freeze();
			for (int iter=0 ; iter<8 ; iter++) {
				DistanceMatrix distanceMatrix = new DistanceMatrix(prot.atoms(), 5.5, 2.0,4); 
				TotalEnergy energy = new TotalEnergy(prot, distanceMatrix, energyCreators, commands);
				Minimizer minimizer = new LBFGS(energy, 0.05, 10000, 1000);
				try {
					minimizer.minimize();
				} catch (Exception e1) {
					throw new RuntimeException(e1);
				}
				energy.resetAtomEnergies();
				AbstractEnergy term = energy.getEnergyTerm(new SoftExcludedVol());
				term.evaluateAtoms();
				int badCounter = 0;
				for (int res=0 ; res<prot.residues().size() ; res++) {
					boolean bad = false;
					for (int atomC=0 ; atomC<prot.residues().residueAt(res).atoms().size() ; atomC++) {
						if (!prot.residues().residueAt(res).atoms().atomAt(atomC).isBackbone &&
								(prot.residues().residueAt(res).atoms().atomAt(atomC).energy()>5.0)) {
							bad = true;
						}
					}
					if (bad) {
						badCounter++;
						double[] tmprot={-999,-999,-999,-999};
						tmprot[0] = 2*Math.PI*(0.5-Math.random());
						tmprot[1] = 2*Math.PI*(0.5-Math.random());
						tmprot[2] = 2*Math.PI*(0.5-Math.random());
						tmprot[3] = 2*Math.PI*(0.5-Math.random());
						ResidueBuilder.build(prot.residues().residueAt(res),prot.residues().residueAt(res).type,tmprot);
					}
				}
				System.out.println("Number of bads: " + badCounter);
			}
			try {
				prot.atoms().print(new MeshiWriter(names[doProt]+".jumbleMin.pdb"));
				System.out.println("Wrote file successfully.");
			} catch (IOException e) {
				throw new RuntimeException("Could not write file.");
			}						
		}
	}

}
