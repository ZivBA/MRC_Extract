package meshi.applications.loopBuilding;

import java.util.Vector;

import meshi.applications.corpus.Corpus;
import meshi.geometry.ResidueBuilder;
import meshi.molecularElements.Atom;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;

public class LoopBuilderUserDefinedFragmentsStochasticWithNCtypeClosure extends LoopBuilderUserDefinedFragmentsStochastic {
	
	protected final double GAP_CUTOFF=8.0; // The building tree is pruned if the gap score is above this value 	
	
	private double closureToleranceSquared = -99999.0;  // We precalculate the square to save computaion time
	private double closureToleranceSquaredIntial=-1; 
	private Atom atom1N,atom1CA,atom1C,atom2N,atom2CA,atom2C; // The closing atoms over the final peptide bond
	int closureRes=-1;
	private double[] xyza = new double[4]; // An auxilary array
	
	// Stats
	private int numberReachClosure = 0;
	private int DISaccepted = 0;
	private int DISrejected = 0;

	public LoopBuilderUserDefinedFragmentsStochasticWithNCtypeClosure(CommandList commands,
			String writePath, Corpus corpus, Protein prot, Protein ref,
			int resStart, int resEnd, double rmsMatchCO, double rmsCutOff,
			Vector<int[]> fragsDescription, double closureTolerance) {
		super(commands, writePath, corpus, prot, ref, resStart, resEnd,
				rmsMatchCO, rmsCutOff, fragsDescription);
		this.closureToleranceSquared = closureTolerance*closureTolerance;
		closureToleranceSquaredIntial = (1.33+closureTolerance)*(1.33+closureTolerance);
		System.out.println("Closure tolerance: " + closureTolerance);
	}

	protected void setClosureAtoms() {
		int lastInd = libs.length - 1;
		if (libManners[lastInd]==-1) // Ends from the N direction
			closureRes = resStart + libEnds[lastInd];
		else   // Ends from the C direction
			closureRes = resStart + libStarts[lastInd]-1;
		atom1C = prot.atoms().findAtomInList("C", closureRes);
		atom2N = prot.atoms().findAtomInList("N", closureRes+1);
		// For the clamp-style closure:
		atom1N = prot.atoms().findAtomInList("N", closureRes);
		atom1CA = prot.atoms().findAtomInList("CA", closureRes);
		atom2CA = prot.atoms().findAtomInList("CA", closureRes+1);
		atom2C = prot.atoms().findAtomInList("C", closureRes+1);

		System.out.println("Closure atoms:\n"+atom1C+"\n"+atom2N);
		System.out.print("Closure residues:\n"+atom1C.residueNumber()+" 9999999999\n"+atom2N.residueNumber()+" 9999999999\n");		
	}
	
	
	protected boolean loopClosingCondition() {
		numberReachClosure++;
		// First: the N and C must be within the tolerance 
		if (((atom1C.x()-atom2N.x())*(atom1C.x()-atom2N.x()) + 
			 (atom1C.y()-atom2N.y())*(atom1C.y()-atom2N.y()) + 
			 (atom1C.z()-atom2N.z())*(atom1C.z()-atom2N.z())) < closureToleranceSquaredIntial) {
			// Second: the extrapolation according to psi of the 2N should be near the real 2N
			ResidueBuilder.getAtom_xyza(xyza, 1.33, 2.037, pp[closureRes][1] , atom1C, atom1CA, atom1N);
			if (((xyza[0]-atom2N.x())*(xyza[0]-atom2N.x()) + 
				 (xyza[1]-atom2N.y())*(xyza[1]-atom2N.y()) + 
				 (xyza[2]-atom2N.z())*(xyza[2]-atom2N.z())) < closureToleranceSquared) {
				// Third: the extrapolation according to phi of the 1C should be near the real 1C
				ResidueBuilder.getAtom_xyza(xyza, 1.33, 2.118, pp[closureRes+1][0] , atom2N, atom2CA, atom2C);
				if (((atom1C.x()-xyza[0])*(atom1C.x()-xyza[0]) + 
					 (atom1C.y()-xyza[1])*(atom1C.y()-xyza[1]) + 
					 (atom1C.z()-xyza[2])*(atom1C.z()-xyza[2])) < closureToleranceSquared) {
					return true;
				}				
			}
		}
		return false;
	}
		
//	protected boolean loopClosingCondition() {
//		numberReachClosure++;
//		// First: the N and C must be within the tolerance 
//		if (((atom1C.x()-atom2N.x())*(atom1C.x()-atom2N.x()) + 
//			 (atom1C.y()-atom2N.y())*(atom1C.y()-atom2N.y()) + 
//			 (atom1C.z()-atom2N.z())*(atom1C.z()-atom2N.z())) < closureToleranceSquared) {
//			// Second: the extrapolation according to psi of the 2N should be near the real 2N
//			ResidueBuilder.getAtom_xyza(xyza, 1.33, 2.037, pp[closureRes][1] , atom1C, atom1CA, atom1N);
//			if (((xyza[0]-atom2N.x())*(xyza[0]-atom2N.x()) + 
//				 (xyza[1]-atom2N.y())*(xyza[1]-atom2N.y()) + 
//				 (xyza[2]-atom2N.z())*(xyza[2]-atom2N.z())) < closureToleranceSquared) {
//				// Third: the extrapolation according to phi of the 1C should be near the real 1C
//				ResidueBuilder.getAtom_xyza(xyza, 1.33, 2.118, pp[closureRes+1][0] , atom2N, atom2CA, atom2C);
//				if (((atom1C.x()-xyza[0])*(atom1C.x()-xyza[0]) + 
//					 (atom1C.y()-xyza[1])*(atom1C.y()-xyza[1]) + 
//					 (atom1C.z()-xyza[2])*(atom1C.z()-xyza[2])) < closureToleranceSquared) {
//					return true;
//				}				
//			}
//		}
//		return false;
//	}

	@Override
	protected boolean continueForNextCall(int libCounter) {
		boolean goodDis = (temporaryDistEnergy(startResForClosingPotential[libCounter],endResForClosingPotential[libCounter])<GAP_CUTOFF);
		if (goodDis)
			DISaccepted++;
		else
			DISrejected++;
		return (goodDis && super.continueForNextCall(libCounter));
	}

	protected void 	printStats() {
	    System.out.println("Number of closure assessed: " + numberReachClosure);	
	    System.out.println("Number of distance acceptances: " + DISaccepted);	
	    System.out.println("Number of distance rejectences: " + DISrejected);	
	    System.out.println("Number of EV acceptances: " + EVaccepted);	
	    System.out.println("Number of EV rejectences: " + EVrejected);	
	    System.out.println("Number of TETHER acceptances: " + TETHERaccepted);	
	    System.out.println("Number of TETHER rejectences: " + TETHERrejected);		    
	}
		
}
