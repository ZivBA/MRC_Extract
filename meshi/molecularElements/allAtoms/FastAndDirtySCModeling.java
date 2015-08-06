package meshi.molecularElements.allAtoms;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.softExcludedVol.SoftExcludedVol;
import meshi.energy.softExcludedVol.SoftExcludedVolCreator;
import meshi.energy.solvate.SolvateCreatorLongHB;
import meshi.energy.solvate.SolvateEnergy;
import meshi.energy.torsionVal.TorsionValEnergy;
import meshi.energy.twoTorsions.TwoTorsionsEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.ResidueBuilder;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionList;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.KeyWords;


public class FastAndDirtySCModeling implements Residues, 
					       AtomTypes , KeyWords{ /**
								      * The implemented
								      * interfaces defines the 
								      * names of atom and residue 
								      * types. 
								      **/
    private DunbrackLib lib;
    private Protein protein; 
    private DistanceMatrix distanceMatrix;
    private TotalEnergy energy;
    //private double min1,min2,min3;
    private int[] seder = {18,19,6,4,9,10,3,8,14,2,13,11,1,12,16,7,17,15}; 
    private double[] w1 = {0.00,0.08,1.78,0.22,2.42,0.00,1.40,0.08,0.50,
			   0.72,0.92,0.32,0.30,0.56,0.30,0.10,0.14,0.04,8.54,7.56};
    private double[] w2 = {0.00,0.02,0.74,0.00,1.04,0.00,1.20,0.10,0.00,
			   0.00,0.46,0.42,0.00,0.00,0.00,0.02,0.04,0.06,1.30,1.22};


    public FastAndDirtySCModeling(CommandList commands) {
	lib = new DunbrackLib(commands, 0.95,20);	
	protein = null;
	distanceMatrix = null;
	energy = null;
    }

    public void model(Protein protein,CommandList commands)  {
	this.protein = protein;
	protein.defrost();       
	protein.atoms().freeze(new AtomList.BackboneFilter());    

	// The energy creators that will be used from now on!
	EnergyCreator[] energyCreators1 = {
	    new SoftExcludedVolCreator(1.0),
	    new SolvateCreatorLongHB(1.0)
	};

	distanceMatrix = new DistanceMatrixWhichUpdateOnlyMovingAtoms(protein.atoms(), 5.5,4); 
	energy = new TotalEnergy(protein, distanceMatrix, energyCreators1, commands);

	int ind;
	int pred;
	int numOfRotamers;
	double[] tmprot;
	double probrot,maxprob;
	double score,minscore;


	double[][] pp = new double[5000][3];
	double[][] targets = new double[5000][4];
	for (int c=0 ; c<targets.length ; c++)
	    targets[c][0] = targets[c][1] = targets[c][2] = targets[c][3] = -999;
	double[] badPoints;
	double[] tmpBadPoints;
	double[][] tmpchis;
	// Reading the phipsi
	phipsi(pp);

	// putting Rot1 into targets
	for (int kk=0 ; kk<pp.length ; kk++) {
	    if ((pp[kk][2]>-900) && (pp[kk][2]!=0) && (pp[kk][2]!=5)) {
		tmprot = lib.getRotamer((int) pp[kk][2] , pp[kk][0] , pp[kk][1] , 0);
		for (int dd=0 ; dd<lib.getChiMax((int) pp[kk][2]) ; dd++)
		    targets[kk][dd] = tmprot[dd];
	    }
	}
	for (int kk=0 ; kk<protein.residues().size() ; kk++) 
	    if (protein.residues().residueAt(kk).type > -1) {
		ResidueBuilder.build(protein.residues().residueAt(kk), 
				     protein.residues().residueAt(kk).type, 
				     targets[kk]);
	    }
   

	// Looping on all the residues.
	for (int cc=0 ; cc<seder.length ; cc++)
	    for (int c=0 ; c<pp.length ; c++) 
		if (pp[c][2]==seder[cc]) {
    
        
		    //*****************************************************************************    
		    // Going over the rotamers
		    numOfRotamers = lib.getRotamerNum((int) pp[c][2] , pp[c][0] , pp[c][1]);
		    badPoints = new double[4*(numOfRotamers)];
		    tmpchis = new double[(numOfRotamers)][4];
    
		    for (ind=0 ; ind<numOfRotamers ; ind++) {
			System.out.println("\nRot " + ind);
			tmprot = lib.getRotamer((int) pp[c][2] , pp[c][0] , pp[c][1] , ind);
			probrot = lib.getRotamerProb((int) pp[c][2] , pp[c][0] , pp[c][1] , ind);
			for (int dd=0 ; dd<lib.getChiMax((int) pp[c][2]) ; dd++)
			    targets[c][dd] = tmprot[dd];

			for (int kk=0 ; kk<protein.residues().size() ; kk++) 
			    if (protein.residues().residueAt(kk).type > -1) {
				ResidueBuilder.build(protein.residues().residueAt(kk), 
						     protein.residues().residueAt(kk).type, 
						     targets[kk]);
			    }
    
			tmpchis[ind][0] = targets[c][0];
			tmpchis[ind][1] = targets[c][1];   
			tmpchis[ind][2] = targets[c][2];
			tmpchis[ind][3] = targets[c][3];   


			// Analyzing the merits of each rotamer                       
			tmpBadPoints = brain(c,probrot); 
			for (int dd=0 ; dd<4 ; dd++) 
			    badPoints[ind*4+dd] = tmpBadPoints[dd];
        
		    } 
		    // END - Going over the rotamers
		    //--------------------------------------------------------------------

		    // Scoring each rotamer
		    minscore = 1e10;
		    pred = -1;
		    maxprob = 0.0;
		    for (ind=0 ; ind<numOfRotamers ; ind++) 
			if (badPoints[ind*4+1]>maxprob)
			    maxprob = badPoints[ind*4+1];
		    for (ind=0 ; ind<numOfRotamers ; ind++) {
			score = badPoints[ind*4] +
			    w1[seder[cc]]*maxprob/(0.01+badPoints[ind*4+1]) +
			    w2[seder[cc]]*badPoints[ind*4+3];
			if (score<minscore) {
			    pred = ind;
			    minscore = score;
			}
		    }
		    if (pred<0)
			throw new RuntimeException("Serious problems " + pred);
           

		    // Returning everything to normal, with the residue in its new rotamer.
		    targets[c][0] = tmpchis[pred][0];
		    targets[c][1] = tmpchis[pred][1];
		    targets[c][2] = tmpchis[pred][2];
		    targets[c][3] = tmpchis[pred][3];
		    for (int resc=0 ; resc<protein.residues().size() ; resc++) 
			if (protein.residues().residueAt(resc).type > -1) {
			    ResidueBuilder.build(protein.residues().residueAt(resc), 
						 protein.residues().residueAt(resc).type, 
						 targets[resc]);
			}    
		} // Going for the next residue.
    }  // of model()



    protected double[] brain(int resnum , double prob) {  

	double[] result = new double[4];
	for (int c=0 ; c<4 ; c++)
	    result[c] = 0;

	offSOL(energy);
	onEV(energy);
	energy.evaluate();
	energy.evaluateAtoms();
	for (int c=0 ; c<protein.atoms().size() ; c++) 
	    if ((protein.atoms().atomAt(c).residueNumber()==resnum) &&
		(protein.atoms().atomAt(c).name().length()>1) && (!protein.atoms().atomAt(c).name().equals("CA")) && 
		(!protein.atoms().atomAt(c).name().equals("CB"))) {
		result[0] += protein.atoms().atomAt(c).energy();
	    }    
      
	result[1] = prob;
    
	onSOL(energy);
	energy.evaluate();
   	AbstractEnergy energyTerm = energy.getEnergyTerm(new SolvateEnergy());
	result[3] = energyTerm.evaluate();


	return result;
    }


    protected void phipsi(double[][] pp) {
	for (int c=0 ; c<pp.length ; c++) {
	    pp[c][0] = pp[c][1] = -60*Math.PI/180.0;
	    pp[c][2]= -999;
	}
	TorsionList phiList = (TorsionList) TorsionList.createTorsionList(protein,
									  distanceMatrix).filter(new TorsionList.FilterPhi(),
												 new TorsionList());
	TorsionList psiList = (TorsionList) TorsionList.createTorsionList(protein,
									  distanceMatrix).filter(new TorsionList.FilterPsi(),
												 new TorsionList());
	for(int i=0 ; i<phiList.size() ; i++) {
	    pp[((Torsion) phiList.elementAt(i)).getTorsionResNum()][0] =
		((Torsion) phiList.elementAt(i)).torsion();
	    pp[((Torsion) phiList.elementAt(i)).getTorsionResNum()][2] = Residue.type(
										      ((Torsion) phiList.elementAt(i)).getTorsionResName());
	}
	for(int i=0 ; i<psiList.size() ; i++) {
	    pp[((Torsion) psiList.elementAt(i)).getTorsionResNum()][1] =
		((Torsion) psiList.elementAt(i)).torsion();
	    pp[((Torsion) psiList.elementAt(i)).getTorsionResNum()][2] = Residue.type(
										      ((Torsion) psiList.elementAt(i)).getTorsionResName());
	}
    }

    
    protected  void onSOL(TotalEnergy te) {
        AbstractEnergy[] energyTerms = te.getEnergyTerms(new SolvateEnergy());
        for (int c=0 ; c<energyTerms.length ; c++) 
	    energyTerms[c].on();
    }

    protected  void offSOL(TotalEnergy te) {
        AbstractEnergy[] energyTerms = te.getEnergyTerms(new SolvateEnergy());
        for (int c=0 ; c<energyTerms.length ; c++) 
	    energyTerms[c].off();
    }

    protected  void onEV(TotalEnergy te) {
        AbstractEnergy[] energyTerms = te.getEnergyTerms(new SoftExcludedVol());
        for (int c=0 ; c<energyTerms.length ; c++) 
	    energyTerms[c].on();
    }

    protected  void offEV(TotalEnergy te) {
        AbstractEnergy[] energyTerms = te.getEnergyTerms(new SoftExcludedVol());
        for (int c=0 ; c<energyTerms.length ; c++) 
	    energyTerms[c].off();
    }

    protected  void on2T(TotalEnergy te) {
        AbstractEnergy[] energyTerms = te.getEnergyTerms(new TwoTorsionsEnergy());
        for (int c=0 ; c<energyTerms.length ; c++) 
	    energyTerms[c].on();
    }

    protected  void off2T(TotalEnergy te) {
        AbstractEnergy[] energyTerms = te.getEnergyTerms(new TwoTorsionsEnergy());
        for (int c=0 ; c<energyTerms.length ; c++) 
	    energyTerms[c].off();
    }

    protected  void onTOR(TotalEnergy te) {
        AbstractEnergy[] energyTerms = te.getEnergyTerms(new TorsionValEnergy());
        for (int c=0 ; c<energyTerms.length ; c++) 
	    energyTerms[c].on();
    }

    protected  void offTOR(TotalEnergy te) {
        AbstractEnergy[] energyTerms = te.getEnergyTerms(new TorsionValEnergy());
        for (int c=0 ; c<energyTerms.length ; c++) 
	    energyTerms[c].off();
    }


} // of class

