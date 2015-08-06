package meshi.IMPmy;

import java.util.Random;

import meshi.molecularElements.AtomList;
import meshi.util.crossLinking.CrosslinkVectorTRiC;

public class TRiC_1_ring extends DomainListWithStructure {

	private static final long serialVersionUID = 1L;

	public TRiC_1_ring(String fileName) {
		super(fileName);
	}

	@Override
	AtomList getAtomicModelFileName(String protName, String domainName) {
		return (new AtomList("One_Ring_Model_OMS.pdb")).chainFilter(protName);
//		return new AtomList("HM_"+protName+".pdb");
	}


	public static void main(String[] args) {
		new Random(2989);
		DomainListWithStructure domList = new TRiC_1_ring("domainsBovine.txt");
		System.out.print("import _surface\n" +
		"import chimera\n" +
		"try:\n" +
		"  import chimera.runCommand\n" +
		"except:\n" +
		"  pass\n" +
		"from VolumePath import markerset as ms\n" +
		"try:\n" +
		"  from VolumePath import Marker_Set, Link\n" +
		"  new_marker_set=Marker_Set\n" +
		"except:\n" +
		"  from VolumePath import volume_path_dialog\n" +
		"  d= volume_path_dialog(True)\n" +
		"  new_marker_set= d.new_marker_set\n" +
		"marker_sets={}\n" +
		"surf_sets={}\n");
		for (Domain dom : domList) {
			System.out.print("if \"\" not in marker_sets:\n" +
					"  s=new_marker_set('')\n" +
					"  marker_sets[\"\"]=s\n" +
					"s= marker_sets[\"\"]\n");
			double red=0.0;
			double green = 0.0;
			double blue = 0.0;
			switch (dom.proteinName().charAt(0)) {
			case 'A':
				 red=1.0;
				 break;
			case 'B':
				 green=1.0;
				 break;
			case 'G':
				 blue=1.0;
				 break;
			case 'D':
				 red=1.0;
				 green=1.0;
				 break;
			case 'E':
				 blue=1.0;
				 green=1.0;
				 break;
			case 'Z':
				 red=102/255.0;
				 green=51/255.0;
				 break;
			case 'H':
				 red=1.0;
				 blue=1.0;
				 break;
			case 'Q':
				 red=75/255.0;
				 green=75/255.0;
				 blue=75/255.0;
				 break;
			}
			System.out.print("mark=s.place_marker(("+dom.center().x()+","+dom.center().y()+","+dom.center().z()+"), ("+red+","+green+","+blue+"), "+2.5*Math.pow(dom.numberOfResidues(),0.333)+")\n");
		}
		System.out.print("for k in surf_sets.keys():\n"+
				"  chimera.openModels.add([surf_sets[k]])\n");
		System.exit(0);
		domList.setBoundaries(1.0);
		domList.setXLvecs(new CrosslinkVectorTRiC("Abersold_Bovine_intra_inter.txt",0),1.0);
		domList.putDomainsInRandomCube(300.0);
		domList.rigidify(1.0);
		domList.addEV(10.0);
		TotalEnergy totalEnergy = new TotalEnergy(domList.disConstList());
//		System.out.println(domList.disConstList());
		SteepestDescent minimizer = new SteepestDescent(totalEnergy, domList);
		minimizer.run(20000000, 0.0001, 200000);
		
		//totalEnergy.evaluate(true);
		//System.out.print(totalEnergy);
//		domList.disConstList().zeroConnectivity();
//		domList.disConstList().upRigidity();
		System.out.println(domList.disConstList());
		
		//System.out.println(domList.disConstList());
	}

//	public static void main(String[] args) {
//		new Random(299);
////		CrosslinkVectorTRiC vec = new CrosslinkVectorTRiC("Abersold_Bovine_intra_inter.txt",0);
////		CrosslinkVectorTRiC vec2 = new CrosslinkVectorTRiC("NK_combined_moderate.txt",1);
////		vec.addVec(vec2);
//		CrosslinkVectorTRiC vec = new CrosslinkVectorTRiC("ArtificialSet_1000.txt",1);
//		DomainListFromSequences domList = new DomainListFromSequences("bovineSeqs.txt");
//		domList.setBoundaries(100.0);
//		domList.setXLvecs(vec.filterOutInterUnit(),1.0);
//		domList.seperateProteins(300.0);
//		domList.addEV(100.0);
//		TotalEnergy totalEnergy = new TotalEnergy(domList.disConstList());
////		System.out.println(domList.disConstList());
//		SteepestDescent minimizer = new SteepestDescent(totalEnergy, domList);
//		minimizer.run(20000000, 0.2, 5000);
//		
//		//totalEnergy.evaluate(true);
//		//System.out.print(totalEnergy);
////		domList.disConstList().zeroConnectivity();
////		domList.disConstList().upRigidity();
//		System.out.println(domList.disConstList());
//		
//		domList.resetConstraintList();
//		domList.setBoundaries(100.0);
//		domList.setXLvecs(vec,1.0);
//		domList.addEV(100.0);
//		totalEnergy = new TotalEnergy(domList.disConstList());
//		minimizer = new SteepestDescent(totalEnergy, domList);
//		minimizer.run(20000000, 0.01, 10000);
//		System.out.println(domList.disConstList());
//	}

	
}



// Code for derivative checks:
//double small = 0.000001;
//totalEnergy.snapshot();
//double e0 = totalEnergy.evaluate(true);
//for (AnchorPosition pos : totalEnergy.positions()) {
//	double grad = pos.Fx();
//	pos.addX(small);
//	double e1 = totalEnergy.evaluate(false);
//	totalEnergy.restoreSnapshot();
//	if ( Math.abs(((e0-e1)/small) - grad) /  Math.abs(((e0-e1)/small) + grad) > 0.001) {
//		System.out.println("X: "+((e0-e1)/small) + "    " + grad + "\n" + pos + "\n\n");
//	}
//	grad = pos.Fy();
//	pos.addY(small);
//	e1 = totalEnergy.evaluate(false);
//	totalEnergy.restoreSnapshot();
//	if ( Math.abs(((e0-e1)/small) - grad) /  Math.abs(((e0-e1)/small) + grad) > 0.001) {
//		System.out.println("Y: "+((e0-e1)/small) + "    " + grad + "\n" + pos + "\n\n");
//	}
//	grad = pos.Fy();
//	pos.addZ(small);
//	e1 = totalEnergy.evaluate(false);
//	totalEnergy.restoreSnapshot();
//	if ( Math.abs(((e0-e1)/small) - grad) /  Math.abs(((e0-e1)/small) + grad) > 0.001) {
//		System.out.println("Z: "+((e0-e1)/small) + "    " + grad + "\n" + pos + "\n\n");
//	}
//	
//}
