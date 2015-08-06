package meshi.energy.rot1Pairwise;

import meshi.energy.NonBondedEnergyElement;
import meshi.energy.ROT1solvation.parameters.WeightRepresentativePolars;
import meshi.energy.ROT1solvation.parameters.WeightsHydrophobicSA;
import meshi.geometry.Distance;
import meshi.molecularElements.Atom;

public  class Rot1PairwiseEnergyElement extends NonBondedEnergyElement implements WeightRepresentativePolars,
WeightsHydrophobicSA {
	// --- Auxilary variables used in the evaluation --------
	private Atom atom1, atom2;
	private double weight;
	private Distance distance;
	private Rot1PairwiseParametersList parametersList;
	private double[] EofR;
	private double maxDis;
	private double dis = -1;
	private double SAfactor = -1;
	// ------------------------------------------------------

	
	public  Rot1PairwiseEnergyElement() {}
	
	public  Rot1PairwiseEnergyElement(Rot1PairwiseParametersList parametersList, 
			int type, double weight){
		this.parametersList = parametersList;
		this.weight = weight;
	}

	protected void setAtoms(){
		throw new RuntimeException("setAtoms() may not be used by Rot1PairwiseEnergyElement for efficiency.");
	}

	public void set(Object obj) {
		distance = (Distance) obj;
		atoms = distance.atoms();
		atom1 = distance.atom1();
		atom2 = distance.atom2();
		Rot1PairwiseParameters parameters = (Rot1PairwiseParameters) parametersList.parameters(distance);
		EofR = parameters.EofR;
		maxDis = parameters.maxDis;
	}        

	public double evaluate() {
		dis = distance.distance();
		if (dis>maxDis)
			return 0.0;
		if ((atom1.residueNumber() > (atom2.residueNumber()-2)) &&
				(atom1.residueNumber() < (atom2.residueNumber()+2)))
			return 0.0;
				
		/* Perhaps we need to evaluate derivatives:
		if (! atom1.frozen()) {
			atom1.addToFx(-1*dEdX*tmpW); // force = -derivative   
			atom1.addToFy(-1*dEdY*tmpW);
			atom1.addToFz(-1*dEdZ*tmpW);
		}
		if (! atom2.frozen()) {
			atom2.addToFx(dEdX*tmpW);
			atom2.addToFy(dEdY*tmpW);
			atom2.addToFz(dEdZ*tmpW);
		}
		*/
//		if (((atom1.residueNumber()==87) && atom1.name().equals("CB")) ||
//				((atom2.residueNumber()==87) && atom2.name().equals("CB")))
//			System.out.println(this+""+dis + " " +weight*EofR[(int) Math.round(dis)]);

		
		SAfactor = -999;
		if (atom1.name().equals("CB") && atom2.name().equals("CB")) {
			if ((weightRepresentativePolars[atom1.type]>0.1) && (weightRepresentativePolars[atom2.type]>0.1)) { // Two philic
				SAfactor = 0.000000;
			}
			else if ((weightRepresentativePolars[atom1.type]<0.1) && (weightRepresentativePolars[atom2.type]>0.1)) { // 1 is phobic 2 is philic
				SAfactor = 0.00000*weightsHydrophobicSA[atom1.type];
			}
			else if ((weightRepresentativePolars[atom1.type]>0.1) && (weightRepresentativePolars[atom2.type]<0.1)) { // 1 is philic 2 is phobic
				SAfactor = 0.00000*weightsHydrophobicSA[atom2.type];
			}
			else  { // Two phobic
				SAfactor = Math.sqrt(weightsHydrophobicSA[atom1.type]*weightsHydrophobicSA[atom2.type]);
			}
		}
		if (atom1.name().equals("CB") && atom2.name().equals("CA")) {
			if (weightRepresentativePolars[atom1.type]>0.1) { // 1 is philic
				SAfactor = 0.00000; //0.7;  
			}
			else  { // 1 is phobic
				SAfactor = Math.sqrt(weightsHydrophobicSA[atom1.type]*0.7); 
			}
		}
		if (atom1.name().equals("CA") && atom2.name().equals("CB")) {
			if (weightRepresentativePolars[atom2.type]>0.1) { // 2 is philic
				SAfactor = 0.00000; //0.7;  
			}
			else  { // 2 is phobic
				SAfactor = Math.sqrt(weightsHydrophobicSA[atom2.type]*0.7); 
			}
		}
		if (atom1.name().equals("CA") && atom2.name().equals("CA")) {
			SAfactor = 0.0;
		}
		
		if (SAfactor<-500)
			throw new RuntimeException("Impossible option:\n" + atom1 + "\n" + atom2+"\n");

		int ind = (int) Math.floor(dis);
		if (dis<maxDis)
			return weight*SAfactor*((1-(dis-ind))*EofR[ind] + (dis-ind)*EofR[ind+1]);
		else
			return 0.0;
	}

	public String toString() {
		if ((atom1 == null) & (atom2 == null)) return "Rot1PairwiseEnergyElement - atoms not set yet";
		if ((atom1 == null) | (atom2 == null)) throw new RuntimeException("This is weird\n"+
				"atom1 = "+atom1+"\n"+
				"atom2 = "+atom2); 
		return ("Rot1PairwiseEnergyElement atoms: \n"+atom1+"\n"+atom2+"\n");
	}

}
