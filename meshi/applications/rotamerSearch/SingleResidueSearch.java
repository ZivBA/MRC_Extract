package meshi.applications.rotamerSearch;

import meshi.geometry.ResidueBuilder;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Residue;

public class SingleResidueSearch {
	
	public final int rangeAroundRotamer = 1; // in multiplications of the resolution
	public final double resolutionAroundRotamer = (Math.PI / 180.0) * 15.0;	
//	public final int rangeAroundRotamer = 3; // in multiplications of the resolution
//	public final double resolutionAroundRotamer = (Math.PI / 180.0) * 8.0;	
	private Residue res;
	private double[][] conformations;
	private double[] probabilities;
	private int pointer;
	private int best;
	private double[][] initialCoordinates;
	
	public SingleResidueSearch(Residue res , double phi, double psi, DunbrackLib lib) {
		if ((res.type<1) | (res.type>19) | (res.type==5)) {
			throw new RuntimeException("Unknown residue type or residue without rotamers.");
		}
		this.res = res;
		pointer = -1;
		best = -1;
		int rangeChi1 = rangeAroundRotamer;
		int rangeChi2 = rangeAroundRotamer;
		int rangeChi3 = rangeAroundRotamer;
		int rangeChi4 = rangeAroundRotamer;
		if (lib.getChiMax(res.type)<4) {
			rangeChi4 = 0;
		}
		if (lib.getChiMax(res.type)<3) {
			rangeChi3 = 0;
		}
		if (lib.getChiMax(res.type)<2) {
			rangeChi2 = 0;
		}
		if (lib.getChiMax(res.type)==3) {
			rangeChi1 = 0;
		}
		if (lib.getChiMax(res.type)==4) {
			rangeChi1 = 0;
			rangeChi2 = 0;
		}
		int possibilities = (2*rangeChi1+1)*(2*rangeChi2+1)*(2*rangeChi3+1)*(2*rangeChi4+1)*lib.getRotamerNum(res.type, phi, psi);
		conformations = new double[possibilities][4];
		probabilities = new double[possibilities];
		int counter = 0;
		for (int rotNum = 0 ; rotNum < lib.getRotamerNum(res.type, phi, psi) ; rotNum++) {
			double[] tmpRot = lib.getRotamer(res.type, phi, psi, rotNum);
			double rotProb = lib.getRotamerProb(res.type, phi, psi, rotNum);
			for (int chi1=-rangeChi1 ; chi1<=rangeChi1 ; chi1++) {
				for (int chi2=-rangeChi2 ; chi2<=rangeChi2 ; chi2++) {
					for (int chi3=-rangeChi3 ; chi3<=rangeChi3 ; chi3++) {
						for (int chi4=-rangeChi4 ; chi4<=rangeChi4 ; chi4++) {
							if (tmpRot.length>3) {
								conformations[counter][3] = tmpRot[3]+chi4*resolutionAroundRotamer;
							} else {
								conformations[counter][3] = 0.0;
							}
							if (tmpRot.length>2) {
								conformations[counter][2] = tmpRot[2]+chi3*resolutionAroundRotamer;
							} else {
								conformations[counter][2] = 0.0;
							}							
							if (tmpRot.length>1) {
								conformations[counter][1] = tmpRot[1]+chi2*resolutionAroundRotamer;
							} else {
								conformations[counter][1] = 0.0;
							}
							conformations[counter][0] = tmpRot[0]+chi1*resolutionAroundRotamer;
							probabilities[counter] = rotProb;
							counter++;							
						}
					}
				}
			}
		}
		// Saving the initial coordinates
		initialCoordinates = new double[res.atoms().size()][3];
		for (int c=0 ; c<res.atoms().size() ; c++) {
			initialCoordinates[c][0] = res.atoms().atomAt(c).x();
			initialCoordinates[c][1] = res.atoms().atomAt(c).y();
			initialCoordinates[c][2] = res.atoms().atomAt(c).z();
		}
	}
	
	/**
	 * Puts first possibility on. 
	 **/
	public void intializeSearch() {
		pointer = 0;
		build(pointer);		
	}
	
	/**
	 * Next option. 
	 **/
	public void nextOption() {
		if (pointer!=conformations.length) {
			pointer++;
		}
		if (pointer!=conformations.length) {
			build(pointer);
		}
	}

	/**
	 * Is valid rotamer. 
	 **/
	public boolean isValid() {
		if ((pointer==-1) | (pointer==conformations.length))
			return false;
		else
			return true;
	}

	/**
	 * Total number of rotamers. 
	 **/
	public int howMany() {
		return conformations.length;
	}
	
	/**
	 * Probability of current rot. 
	 **/
	public double getProb() {
		return probabilities[pointer];
	}
	
	
	/**
	 * The best rotamer so far
	 **/
	public void setBest() {
		best = pointer;
	}
	
	public int pointer() {
		return pointer;
	}

	public int resNumber() {
		return res.number;
	}

	public String resName() {
		return res.name;
	}

	
	/**
	 * Put best
	 **/
	public void buildBest() {
		if (best!=-1) {
			build(best);
		}
		else {
			throw new RuntimeException("Best was not set yet");
		}		
	}
	
	public void restoreInitialCoordinates() {
		for (int c=0 ; c<res.atoms().size() ; c++) {
			res.atoms().atomAt(c).setXYZ(initialCoordinates[c][0] , initialCoordinates[c][1] , initialCoordinates[c][2]);
		}
	}


	private void build(int conformationInd) {
		ResidueBuilder.build(res,res.type,conformations[conformationInd]);
	}
	
	public String toString() {
		return "Pointer: " + pointer() + " ---  " + res;
	}
	
}
