package meshi.energy.softExcludedVol;
import meshi.energy.NonBondedEnergyElement;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;


public  class SoftExcludedVolEnergyElement extends NonBondedEnergyElement {
	protected DistanceMatrix distanceMatrix;
	protected Atom atom1, atom2;
	protected Distance distance;
	int type;
	protected double sigma,C;
	protected boolean frozen;
	protected double dEdD;
	protected double dEdX;
	protected double dEdY;
	protected double dEdZ;
	protected double energy;
	protected double weight;
	protected double rMax;
	protected SoftExcludedVolParametersList parametersList;
	protected int fragStart=100000000,fragEnd=1000000000;

	public  SoftExcludedVolEnergyElement() {}
	public  SoftExcludedVolEnergyElement(SoftExcludedVolParametersList parametersList, DistanceMatrix distanceMatrix,
			int type, double weight) {
		this.parametersList = parametersList;
		this.type = type;
		this.weight = weight;
		this.distanceMatrix = distanceMatrix;
		this.rMax = distanceMatrix.rMax();
		//----------------------
	}
	protected void setAtoms(){
		throw new RuntimeException("setAtoms() may not be used by SoftExcludedVolEnergyElement for "+
		"efficiency.");
	}

	public void setFragLimits(int start , int end) {
		fragStart = start;
		fragEnd = end;
	}


	public void set(Object obj) {
		distance = (Distance) obj;
		atoms = distance.atoms();
		atom1 = distance.atom1();
		atom2 = distance.atom2();
		SoftExcludedVolParameters parameters = (SoftExcludedVolParameters) parametersList.parameters(distance);
		sigma = parameters.sigma;
		//if (!atom1.isBackbone || !atom2.isBackbone)
		//	sigma = sigma*0.8;
		C = parameters.C;
		if (sigma>rMax) 
			throw new RuntimeException("SoftExcluded Volume: sigma="+sigma+" and it is larger " +
					"than rMax="+rMax);
	}


	public double evaluate() {
		double dEdD;
		double DMinusSig;
		double D3;
		double localSigma=sigma;
		double dis = -1;
		dis = distance.distance();
		if (dis>localSigma) {
			energy = 0;
		}
		else {
			DMinusSig = (dis-localSigma);
			if (type==0) {
				if ((atom1.name().length()==1) || (atom2.name().length()==1) || atom1.name().equals("CA") 
						|| atom2.name().equals("CA") || atom1.name().equals("CB") || atom2.name().equals("CB")) {
					D3 = 1.0/(0.2*0.2*0.2*0.2)*DMinusSig*DMinusSig*DMinusSig*weight;
					energy = D3*DMinusSig;
					dEdD = 4*D3;	    	
				}
				else {
					D3 = C*DMinusSig*weight;
					energy = D3*DMinusSig;
					dEdD = 2*D3;	    	
				}
				if (! atom1.frozen()) {
					atom1.addToFx(-1*dEdD*distance.dDistanceDx()); // force = -derivative   
					atom1.addToFy(-1*dEdD*distance.dDistanceDy());
					atom1.addToFz(-1*dEdD*distance.dDistanceDz());
				}
				if (! atom2.frozen()) {
					atom2.addToFx(dEdD*distance.dDistanceDx());
					atom2.addToFy(dEdD*distance.dDistanceDy());
					atom2.addToFz(dEdD*distance.dDistanceDz());
				}
			}
			else if (type==1) {
				D3 = DMinusSig*weight;
				energy = D3*DMinusSig;
				dEdD = 2*D3;

				if (! atom1.frozen()) {
					atom1.addToFx(-1*dEdD*distance.dDistanceDx()); // force = -derivative   
					atom1.addToFy(-1*dEdD*distance.dDistanceDy());
					atom1.addToFz(-1*dEdD*distance.dDistanceDz());
				}
				if (! atom2.frozen()) {
					atom2.addToFx(dEdD*distance.dDistanceDx());
					atom2.addToFy(dEdD*distance.dDistanceDy());
					atom2.addToFz(dEdD*distance.dDistanceDz());
				}
			}
			else if (type==2) {
				if (((atom1.name().length()==1)  || atom1.name().equals("CA") || atom1.name().equals("CB")) &&
						(atom2.name().equals("CB") || atom2.name().equals("CA") || (atom2.name().length()==1))) {
					D3 = 10*DMinusSig*weight;
					energy = D3*DMinusSig;
					dEdD = 2*D3;
				}
				else {
					D3 = DMinusSig*weight;
					energy = D3*DMinusSig;
					dEdD = 2*D3;
				}
				if (! atom1.frozen()) {
					atom1.addToFx(-1*dEdD*distance.dDistanceDx()); // force = -derivative   
					atom1.addToFy(-1*dEdD*distance.dDistanceDy());
					atom1.addToFz(-1*dEdD*distance.dDistanceDz());
				}
				if (! atom2.frozen()) {
					atom2.addToFx(dEdD*distance.dDistanceDx());
					atom2.addToFy(dEdD*distance.dDistanceDy());
					atom2.addToFz(dEdD*distance.dDistanceDz());
				}
			}
			else if (type==3) {
				D3 = 1.0*DMinusSig*weight;
				energy = D3*DMinusSig;
				dEdD = 2*D3;	    	
				if (! atom1.frozen()) {
					atom1.addToFx(-1*dEdD*distance.dDistanceDx()); // force = -derivative   
					atom1.addToFy(-1*dEdD*distance.dDistanceDy());
					atom1.addToFz(-1*dEdD*distance.dDistanceDz());
				}
				if (! atom2.frozen()) {
					atom2.addToFx(dEdD*distance.dDistanceDx());
					atom2.addToFy(dEdD*distance.dDistanceDy());
					atom2.addToFz(dEdD*distance.dDistanceDz());
				}
			}
			else if (type==4) {  // A clash between backbone and non-backbone
				if ((((atom1.name().length()==1) || atom1.name().equals("CA") || atom1.name().equals("CB")) &&
						!((atom2.name().length()==1) || atom2.name().equals("CA") || atom2.name().equals("CB"))) ||
						!(((atom1.name().length()==1) || atom1.name().equals("CA") || atom1.name().equals("CB")) &&
								((atom2.name().length()==1) || atom2.name().equals("CA") || atom2.name().equals("CB"))))   {
					D3 = DMinusSig*weight;
					energy = D3*DMinusSig;
					dEdD = 2*D3;	    	
				}
				else {
					energy = 0;
					dEdD = 0;	    	
				}
				if (! atom1.frozen()) {
					atom1.addToFx(-1*dEdD*distance.dDistanceDx()); // force = -derivative   
					atom1.addToFy(-1*dEdD*distance.dDistanceDy());
					atom1.addToFz(-1*dEdD*distance.dDistanceDz());
				}
				if (! atom2.frozen()) {
					atom2.addToFx(dEdD*distance.dDistanceDx());
					atom2.addToFy(dEdD*distance.dDistanceDy());
					atom2.addToFz(dEdD*distance.dDistanceDz());
				}
			}
			else if (type==9) {  // At least one of the atoms is from the loop (fragment). 
				if (atom1.isBackbone && atom2.isBackbone  &&
						( ((atom1.residueNumber()<=fragEnd)&&(atom1.residueNumber()>=fragStart)) ||
								((atom2.residueNumber()<=fragEnd)&&(atom2.residueNumber()>=fragStart)) ) ) {
					energy =  DMinusSig*DMinusSig*weight;
					dEdD = 2*DMinusSig*weight;
				}
				else {
					energy = 0;
					dEdD = 0;
				}
				if (! atom1.frozen()) {
					atom1.addToFx(-1*dEdD*distance.dDistanceDx()); // force = -derivative
					atom1.addToFy(-1*dEdD*distance.dDistanceDy());
					atom1.addToFz(-1*dEdD*distance.dDistanceDz());
				}
				if (! atom2.frozen()) {
					atom2.addToFx(dEdD*distance.dDistanceDx());
					atom2.addToFy(dEdD*distance.dDistanceDy());
					atom2.addToFz(dEdD*distance.dDistanceDz());
				}
			}
			else if (type==10) {  // One of the atoms is from the loop (fragment) and the other is NOT.
				if (atom1.isBackbone && atom2.isBackbone  &&
						( ((atom1.residueNumber()<=fragEnd)&&(atom1.residueNumber()>=fragStart)&&
								((atom2.residueNumber()>fragEnd)||(atom2.residueNumber()<fragStart)))
								|| // Of the switch between atom 1 and 2
								((atom2.residueNumber()<=fragEnd)&&(atom2.residueNumber()>=fragStart)&&
										((atom1.residueNumber()>fragEnd)||(atom1.residueNumber()<fragStart))) ) ) {
					energy =  DMinusSig*DMinusSig*weight;
					dEdD = 2*DMinusSig*weight;
				}
				else {
					energy = 0;
					dEdD = 0;
				}
				if (! atom1.frozen()) {
					atom1.addToFx(-1*dEdD*distance.dDistanceDx()); // force = -derivative
					atom1.addToFy(-1*dEdD*distance.dDistanceDy());
					atom1.addToFz(-1*dEdD*distance.dDistanceDz());
				}
				if (! atom2.frozen()) {
					atom2.addToFx(dEdD*distance.dDistanceDx());
					atom2.addToFy(dEdD*distance.dDistanceDy());
					atom2.addToFz(dEdD*distance.dDistanceDz());
				}
			}
			else if (type==11) { // This type is for soft backbone EV, that is HARD on HB forming atoms (N,O)
				if ((atom1.name().equals("O") && atom2.name().equals("N")) ||
						(atom1.name().equals("N") && atom2.name().equals("O"))) {
					D3 = 10.0*DMinusSig*weight;
					energy = D3*DMinusSig;
					dEdD = 2*D3;	    	
				}
				else if (((atom1.name().length()==1)  || atom1.name().equals("CA") || atom1.name().equals("CB")) &&
						(atom2.name().equals("CB") || atom2.name().equals("CA") || (atom2.name().length()==1))) {
					D3 = 1.0*DMinusSig*weight;
					energy = D3*DMinusSig;
					dEdD = 2*D3;	    	
				}
				else {
					energy = 0.0;
					dEdD = 0.0;	    	
				}
				if (! atom1.frozen()) {
					atom1.addToFx(-1*dEdD*distance.dDistanceDx()); // force = -derivative   
					atom1.addToFy(-1*dEdD*distance.dDistanceDy());
					atom1.addToFz(-1*dEdD*distance.dDistanceDz());
				}
				if (! atom2.frozen()) {
					atom2.addToFx(dEdD*distance.dDistanceDx());
					atom2.addToFy(dEdD*distance.dDistanceDy());
					atom2.addToFz(dEdD*distance.dDistanceDz());
				}
			}
			else if (type==12) {       //
				if (((atom1.name().length()==1)  || atom1.name().equals("CA") || atom1.name().equals("CB")) &&
						(atom2.name().equals("CB") || atom2.name().equals("CA") || (atom2.name().length()==1))) {
					D3 = DMinusSig*weight;
					energy = D3*DMinusSig;
					dEdD = 2*D3;
				}
				else {
					energy = 0;
					dEdD = 0;
				}
				if (! atom1.frozen()) {
					atom1.addToFx(-1*dEdD*distance.dDistanceDx()); // force = -derivative   
					atom1.addToFy(-1*dEdD*distance.dDistanceDy());
					atom1.addToFz(-1*dEdD*distance.dDistanceDz());
				}
				if (! atom2.frozen()) {
					atom2.addToFx(dEdD*distance.dDistanceDx());
					atom2.addToFy(dEdD*distance.dDistanceDy());
					atom2.addToFz(dEdD*distance.dDistanceDz());
				}
			}
			else
				throw new RuntimeException("An unknown excluded volume type");
		}	
		return energy;
	}


	public String toString() {
		if ((atom1 == null) & (atom2 == null)) return "SoftExcludedVolumeEnergyElement - atoms not set yet";
		if ((atom1 == null) | (atom2 == null)) throw new RuntimeException("This is weird\n"+
				"atom1 = "+atom1+"\n"+
				"atom2 = "+atom2); 
		double dis = distance.distance();
		return ("SoftExcludedVolumeEnergyElement sigma = "+sigma+" Distance = "+
				dFormatSrt.f(dis)+" rMax = "+rMax+"\n"+atom1.verbose(1)+"\n"+atom2.verbose(1));
	}
}
