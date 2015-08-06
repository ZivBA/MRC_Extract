package meshi.applications.TriC.yeast;

import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;

public class AssembleParticleFromSCWRL {

	public static void main(String[] args) {

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter("refine_16-half.ENCAD.SCWRL.pdb"));

			String chain = "A";
			AtomList allAtoms = new AtomList("do_"+chain+"_scwrl.pdb"); Atom.resetNumberOfAtoms();
			AtomList justChain = allAtoms.chainFilter(chain);
			for (int c=0 ; c<justChain.size() ; c++) {
				bw.write(justChain.atomAt(c).toString() + "\n");
			}
			bw.write("TER\n");
			chain = "B";
			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
			justChain = allAtoms.chainFilter(chain);
			for (int c=0 ; c<justChain.size() ; c++) {
				bw.write(justChain.atomAt(c).toString() + "\n");
			}
			bw.write("TER\n");
			chain = "G";
			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
			justChain = allAtoms.chainFilter(chain);
			for (int c=0 ; c<justChain.size() ; c++) {
				bw.write(justChain.atomAt(c).toString() + "\n");
			}
			bw.write("TER\n");
			chain = "D";
			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
			justChain = allAtoms.chainFilter(chain);
			for (int c=0 ; c<justChain.size() ; c++) {
				bw.write(justChain.atomAt(c).toString() + "\n");
			}
			bw.write("TER\n");
			chain = "E";
			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
			justChain = allAtoms.chainFilter(chain);
			for (int c=0 ; c<justChain.size() ; c++) {
				bw.write(justChain.atomAt(c).toString() + "\n");
			}
			bw.write("TER\n");
			chain = "H";
			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
			justChain = allAtoms.chainFilter(chain);
			for (int c=0 ; c<justChain.size() ; c++) {
				bw.write(justChain.atomAt(c).toString() + "\n");
			}
			bw.write("TER\n");
			chain = "Q";
			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
			justChain = allAtoms.chainFilter(chain);
			for (int c=0 ; c<justChain.size() ; c++) {
				bw.write(justChain.atomAt(c).toString() + "\n");
			}
			bw.write("TER\n");
			chain = "Z";
			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
			justChain = allAtoms.chainFilter(chain);
			for (int c=0 ; c<justChain.size() ; c++) {
				bw.write(justChain.atomAt(c).toString() + "\n");
			}			
//			String chain = "A";
//			AtomList allAtoms = new AtomList("do_"+chain+"_scwrl.pdb"); Atom.resetNumberOfAtoms();
//			AtomList justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");
//			chain = "B";
//			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
//			justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");
//			chain = "C";
//			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
//			justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");
//			chain = "D";
//			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
//			justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");
//			chain = "E";
//			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
//			justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");
//			chain = "F";
//			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
//			justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");
//			chain = "G";
//			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
//			justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");
//			chain = "H";
//			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
//			justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");
//			chain = "I";
//			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
//			justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");
//			chain = "J";
//			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
//			justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");
//			chain = "K";
//			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
//			justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");
//			chain = "L";
//			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
//			justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");
//			chain = "M";
//			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
//			justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");
//			chain = "N";
//			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
//			justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");
//			chain = "O";
//			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
//			justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
//			bw.write("TER\n");
//			chain = "P";
//			allAtoms = new AtomList("do_"+chain+"_scwrl.pdb");  Atom.resetNumberOfAtoms();
//			justChain = allAtoms.chainFilter(chain);
//			for (int c=0 ; c<justChain.size() ; c++) {
//				bw.write(justChain.atomAt(c).toString() + "\n");
//			}
			bw.write("TER\n");
			bw.write("END\n");
			bw.close();
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}
	}

	
}
