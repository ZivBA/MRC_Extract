package meshi.energy.caAttraction;
import meshi.energy.CooperativeEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;

public class CaAttractionEnergy extends CooperativeEnergyTerm {
	
/**
 * A inverse square attraction between Ca's that are farther than "APART" on the chain. The attraction is
 * quenced at "rMax" with a margin of "ALPHA"
 **/

    private final int APART = 8;
    private final double rMax = 9.0;
    private final double ALPHA = 0.5;
    private final double rMaxSquared = rMax*rMax;
	private Atom[][] caPairs = null;

    public CaAttractionEnergy() {}
    
    public CaAttractionEnergy(AtomList atomList, 
                    DistanceMatrix dm,
				    double weight) {
	super(toArray(),atomList, dm, null , weight);
	comment = "CaAttraction";
	int counter=0;
	for(int c=0; c<atomList.size() ; c++) {
		for(int d=0; d<atomList.size() ; d++) {
			if ((Math.abs(atomList.atomAt(c).residueNumber()-atomList.atomAt(d).residueNumber()) > APART) &&
			    (!atomList.atomAt(c).frozen() || !atomList.atomAt(d).frozen()))
				counter++;
		}
	}
	caPairs = new Atom[counter][2];
	counter=0;
	for(int c=0; c<atomList.size() ; c++) {
		for(int d=0; d<atomList.size() ; d++) {
			if ((Math.abs(atomList.atomAt(c).residueNumber()-atomList.atomAt(d).residueNumber()) > APART) &&
			    (!atomList.atomAt(c).frozen() || !atomList.atomAt(d).frozen())) {
				caPairs[counter][0] = atomList.atomAt(c);
				caPairs[counter][1] = atomList.atomAt(d);
				counter++;
			}
		}
	}	
    }
 
    public void evaluateAtoms() {
     	evaluate(true);
    }
    
    public double evaluate() {
    	return evaluate(false);
    }
    
    public double evaluate(boolean evaluateAtoms) {
	if (!on) return 0.0;
 	double dEfulldD,Efull,dEdD,energy = 0;
 	double dis,dis2,invdis,dx,dy,dz;
    double rMaxMinusDis,rMaxMinusDisSquare,rMaxMinusDisSquarePlusAlpha,
           rMaxMinusDisTimesAlpha,rMaxMinusDisSquarePlusAlphaSquare,contact,dCdD; 
	
 	for (int c=0 ; c<caPairs.length ; c++) {
 		dx = (caPairs[c][0].x()-caPairs[c][1].x());
 		dy = (caPairs[c][0].y()-caPairs[c][1].y());
 		dz = (caPairs[c][0].z()-caPairs[c][1].z());
 		dis2 = dx*dx + dy*dy + dz*dz;
 	    if (dis2<rMaxSquared) {
 	    	dis = Math.sqrt(dis2);
 	    	invdis = 1/dis;
 	    	rMaxMinusDis = rMax - dis;
            Efull = -invdis;                      // here is the fuctional form of the attraction
            dEfulldD = invdis * invdis;           // ... and its derivative 
            rMaxMinusDisSquare = rMaxMinusDis*rMaxMinusDis;
            rMaxMinusDisSquarePlusAlpha = rMaxMinusDisSquare+ALPHA;
            rMaxMinusDisTimesAlpha = rMaxMinusDis*ALPHA;
            rMaxMinusDisSquarePlusAlphaSquare = rMaxMinusDisSquarePlusAlpha*rMaxMinusDisSquarePlusAlpha;
            contact = rMaxMinusDisSquare/rMaxMinusDisSquarePlusAlpha;
            dCdD = -2*rMaxMinusDisTimesAlpha/rMaxMinusDisSquarePlusAlphaSquare; //The minus comes from rMaxMinusDis
            energy += weight*Efull*contact;
 	    	if (evaluateAtoms) {
 	    		caPairs[c][0].addEnergy(weight*Efull*contact/2); 	    	 
 	    		caPairs[c][1].addEnergy(weight*Efull*contact/2);
 	    	}
            dEdD = weight*(dEfulldD*contact + dCdD*Efull);             

			if (! caPairs[c][0].frozen()) {
				caPairs[c][0].addToFx(-1*dEdD*dx*invdis); // force = -derivative   
				caPairs[c][0].addToFy(-1*dEdD*dy*invdis);
				caPairs[c][0].addToFz(-1*dEdD*dz*invdis);
			}
			if (! caPairs[c][1].frozen()) {
				caPairs[c][1].addToFx(dEdD*dx*invdis);
				caPairs[c][1].addToFy(dEdD*dy*invdis);
				caPairs[c][1].addToFz(dEdD*dz*invdis);
			}
 	    }  	
 	}
 	
    return energy;
    }
}
