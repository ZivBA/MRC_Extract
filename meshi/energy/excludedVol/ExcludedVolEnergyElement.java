package meshi.energy.excludedVol;
import meshi.energy.NonBondedEnergyElement;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;


public  class ExcludedVolEnergyElement extends NonBondedEnergyElement {
    protected DistanceMatrix distanceMatrix;
    protected Atom atom1, atom2;
    protected Distance distance;
    protected double sigma,C,Rfac;
    protected boolean frozen;
    protected double dEdD;
    protected double dEdX;
    protected double dEdY;
    protected double dEdZ;
    protected double energy;
    protected double weight;
    protected double rMax;
    protected ExcludedVolParametersList parametersList;

    public  ExcludedVolEnergyElement() {}
    public  ExcludedVolEnergyElement(ExcludedVolParametersList parametersList, DistanceMatrix distanceMatrix,
               double Rfac, double weight) {
	this.parametersList = parametersList;
	this.weight = weight;
	this.distanceMatrix = distanceMatrix;
	this.Rfac = Rfac;
	this.rMax = distanceMatrix.rMax();
        //----------------------
    }
    protected void setAtoms(){
	throw new RuntimeException("setAtoms() may not be used by LennardJonesEnergyElement for "+
				   "efficiency.");
    }
    
	
    public void set(Object obj) {
	distance = (Distance) obj;
	atoms = distance.atoms();
	atom1 = distance.atom1();
	atom2 = distance.atom2();
	ExcludedVolParameters parameters = (ExcludedVolParameters) parametersList.parameters(distance);
    sigma = parameters.sigma;
    C = parameters.C;
    if (sigma>rMax) 
       throw new RuntimeException("Excluded Volume: sigma="+sigma+" and it is larger " +
                                  "than rMax="+rMax);
    }

        
    public double evaluate() {
	double dEdD;
	double DMinusSig;
	double D3;
	double localSigma=sigma;
	double dis = -1;
	dis = distance.distance();
	if (Rfac<0.99) {
	   if (atom1.residueNumber() == atom2.residueNumber())
	      localSigma = 0;
	   else
	      localSigma = sigma*Rfac;
    }
	if (dis>localSigma) {
	    energy = dEdD = 0;
	}
	else {
	if (((atom1.name().length()==1)  || atom1.name().equals("CA") || atom1.name().equals("CB")) &&
		(atom2.name().equals("CB") || atom2.name().equals("CA") || (atom2.name().length()==1))) {
		DMinusSig = (dis-localSigma);
		D3 = C*DMinusSig*DMinusSig*DMinusSig*weight;
	    energy = D3*(dis-localSigma);
	    dEdD = 4*D3;
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
//if (energy>10.0)
//       System.out.println(atom1.residueNumber+" "+atom2.residueNumber+" "+atom1.type+" "+atom2.type+" "+sigma+" "+dis+" "+A+" "+B+" "+C+" "+ALPHA+" "+sigmaMinusALPHA);
	return energy;
    }
    

    public String toString() {
	if ((atom1 == null) & (atom2 == null)) return "ExcludedVolumeEnergyElement - atoms not set yet";
	if ((atom1 == null) | (atom2 == null)) throw new RuntimeException("This is weird\n"+
									  "atom1 = "+atom1+"\n"+
									  "atom2 = "+atom2); 
	double dis = distance.distance();
	return ("ExcludedVolumeEnergyElement sigma = "+sigma+" Distance = "+
		dFormatSrt.f(dis)+" rMax = "+rMax+"\n"+atom1.verbose(1)+"\n"+atom2.verbose(1));
    }
}
