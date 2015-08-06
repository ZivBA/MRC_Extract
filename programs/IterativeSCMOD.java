package programs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Iterator;
import java.util.StringTokenizer;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.bond.BondCreator;
import meshi.energy.simpleHPterm.SimpleHP;
import meshi.energy.simpleHPterm.SimpleHPCreator;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBondEnergy;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_LowAccuracy_Creator;
import meshi.energy.softExcludedVol.SoftExcludedVolParameters;
import meshi.energy.softExcludedVol.SoftExcludedVolParametersList;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.ResidueBuilder;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.Minimizer;
import meshi.parameters.AtomTypes;
import meshi.parameters.MeshiPotential;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.MeshiList;
import meshi.util.MeshiProgram;
import meshi.util.dssp.DSSP;
import meshi.util.rotamericTools.RotamericTools;


public class IterativeSCMOD extends MeshiProgram implements Residues, 
AtomTypes , KeyWords , MeshiPotential { /**
 * The implemented
 * interfaces defines the 
 * names of atom and residue 
 * types. 
 **/
	private static String commandsFileName = null;
	private static String modelFileName = null;  
	private static String dsspFileName = null;  
	private static String weightsFileName = null;  
	private static String weightsFileName1 = null;  
	private static String writeDir = null;  
	private static Protein fullprot,protein,writeProt; 
	private static int TH1 , TH2;
	private static double A,B,C,D,E;
	private static DecimalFormat fmt = new DecimalFormat("0.#");    
	private static DecimalFormat fmt1 = new DecimalFormat("0.##");    
	private static DecimalFormat fmt2 = new DecimalFormat("0.###");    
	private static DecimalFormat fmt3 = new DecimalFormat("0.####");    
	private static DistanceMatrix distanceMatrix;
	private static 	TotalEnergy energy;
	private static 	SimpleHydrogenBondEnergy hbTerm;
	private static 	SimpleHP solvateTerm;
	private static 	Minimizer minimizer;
	private static double min1,min2,min3;
	private static int[] seder = {18,19,4,9,10,12,17,7,6,3,8,14,2,13,11,1,16,15};


	public static void main(String[] args)  throws Exception {
		CommandList commands; 
		AtomList atoms;
		AbstractEnergy energyTerm;
		MeshiList list;
		PrintWriter pw;
		int typ;
		int ind;
		int pred,pre,post;
		double rotJump = 15*Math.PI/180.0;
		double evalLimit = 40*Math.PI/180.0;
		int eachSide;
//		int numOfRotamers,effectiveNumOfRotamers;
		int rotlikenative=-1;
		int rotInd,flipChi,flipInd,jitterInd1,jitterInd2,maxFlipInd,maxJitterInd1,maxJitterInd2;
		double delta;
		double[] tmprot;
		double[] fourTupple = {-999,-999,-999,-999};
		double[] diffvec = {-999,-999,-999,-999};
		double score,minscore,minrms,diff=999,diff1=999,diff1alt=999;

		init(args); 

		commands = new CommandList(commandsFileName);

		SoftExcludedVolParametersList parametersList = 
			new SoftExcludedVolParametersList(commands.firstWord(PARAMETERS_DIRECTORY).secondWord() + "/" + 
					EXCLUDED_VOL_PARAMETERS,E);

		DunbrackLib lib = new DunbrackLib(commands,1.0,90);

		double[][] w1 = readWeight(weightsFileName);
		double[][] w2 = readWeight(weightsFileName1);
		double[][] w;

		// Creating the proteins.	
		fullprot = new Protein(new AtomList(modelFileName),new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
		for (int cc=0 ; cc<fullprot.atoms().size() ; cc++)
			fullprot.atoms().atomAt(cc).setChain("A");
		protein = new Protein(new AtomList(modelFileName),new ResidueExtendedAtoms(ADD_ATOMS));
		for (int cc=0 ; cc<protein.atoms().size() ; cc++)
			protein.atoms().atomAt(cc).setChain("A");
		writeProt = new Protein(new AtomList(modelFileName),new ResidueExtendedAtoms(ADD_ATOMS));
		for (int cc=0 ; cc<protein.atoms().size() ; cc++)
			writeProt.atoms().atomAt(cc).setChain("A");
		DSSP dssp = new DSSP(dsspFileName);					   

		// Creating the Distance Matrix
		protein.defrost();
		double[][] trueAng = RotamericTools.getSidechainTorsions(protein);
		distanceMatrix = new DistanceMatrix(protein.atoms(), 5.5, 2.0 , 4); 
		RotamericTools.putIntoRot1(protein, distanceMatrix , lib);
		protein.atoms().backbone().freeze();
		protein.atoms().filter(new AtomList.CbFilter()).defrost();
		distanceMatrix = new DistanceMatrix(protein.atoms(), 5.5, 2.0 , 4); 

		// Reading the phipsi
		double[][] pp = RotamericTools.phipsi(protein, distanceMatrix);

		double[][] preBest = new double[pp.length][4];
		double[][] preBestMin = new double[pp.length][4];
		double[][] predictBest = new double[pp.length][4];
		double[][] predictBestMin = new double[pp.length][4];
		double[][] postBest = new double[pp.length][4];
		double[][] postBestMin = new double[pp.length][4];
		for (int c=0 ; c<pp.length ; c++) {
			preBest[c][0] = preBest[c][1] = preBest[c][2] = preBest[c][3] = -999;
			preBestMin[c][0] = preBestMin[c][1] = preBestMin[c][2] = preBestMin[c][3] = -999;
			predictBest[c][0] = predictBest[c][1] = predictBest[c][2] = predictBest[c][3] = -999;
			predictBestMin[c][0] = predictBestMin[c][1] = predictBestMin[c][2] = predictBestMin[c][3] = -999;
			postBest[c][0] = postBest[c][1] = postBest[c][2] = postBest[c][3] = -999;
			postBestMin[c][0] = postBestMin[c][1] = postBestMin[c][2] = postBestMin[c][3] = -999;
		}
		double[][] badPoints;
		double[][][] allrot = new double[pp.length][][];
		double[][] probRot = new double[pp.length][];
		double[][] rmsDevRot = new double[pp.length][];    
		double[][] finalChis;
		double[] finalRMS;
		int[] results;


		// ********************************
		// Manipulations on the library
		// ********************************
		for (int c=0 ; c!=pp.length ; c++)  
			if ((pp[c]!=null) && (pp[c][2]>0) && (pp[c][2]!=5)) {
				typ = (int) pp[c][2];

				// Jittering the rotamers only for PHE,HIS,TRP and TYR
				if ((typ==4) || (typ==6) || (typ==18) || (typ==19))
					eachSide = 1;
				else
					eachSide = 0;

				// Fliping HIS ASN and GLN        
				if ((typ==6) || (typ==11)) {
					maxFlipInd = 2;
					flipChi = 1;
				}
				else if (typ==13) {
					maxFlipInd = 2;
					flipChi = 2;
				}
				else {
					maxFlipInd = 1;
					flipChi = 999;
				}    

				maxFlipInd = 1;
				flipChi = 999;

				if ((typ==1) || (typ==15) || (typ==16) || (typ==17)) {
					maxJitterInd1 = 2*eachSide + 1;
					maxJitterInd2 = 1;
				}
				else {
					maxJitterInd1 = 2*eachSide + 1;
					maxJitterInd2 = 2*eachSide + 1;
				}
				int numOfRotamers = lib.getRotamerNum(typ, pp[c][0] , pp[c][1]);
				int effectiveNumOfRotamers = numOfRotamers*maxJitterInd1*maxJitterInd2*maxFlipInd;// + 1;
				System.out.println("\nOriginal rotamer number:" + numOfRotamers + " efective number of rotamers:" + effectiveNumOfRotamers);
				allrot[c] = new double[effectiveNumOfRotamers][lib.getChiMax(typ)];
				probRot[c] = new double[effectiveNumOfRotamers];
				rmsDevRot[c] = new double[effectiveNumOfRotamers];

				ind=0;
				for (rotInd=0 ; rotInd<numOfRotamers ; rotInd++) {
//					System.out.println("Rot " + rotInd);
					tmprot = lib.getRotamer(typ, pp[c][0] , pp[c][1] , rotInd);
					for (flipInd=0 ; flipInd<maxFlipInd ; flipInd++) 
						for (jitterInd1=0 ; jitterInd1<maxJitterInd1 ; jitterInd1++)     
							for (jitterInd2=0 ; jitterInd2<maxJitterInd2 ; jitterInd2++) {
								probRot[c][ind] = lib.getRotamerProb(typ, pp[c][0] , pp[c][1] , rotInd);
								rmsDevRot[c][ind] = 0.0;
								for (int dd=0 ; dd<lib.getChiMax((int) pp[c][2]) ; dd++) {
									delta = 0.0;
									if (dd==0)
										delta += (jitterInd1-eachSide)*rotJump;
									if (dd==1)
										delta += (jitterInd2-eachSide)*rotJump;
									rmsDevRot[c][ind] += delta*delta;
									if (flipChi==dd)
										delta += Math.PI*flipInd;
									allrot[c][ind][dd] = tmprot[dd] + delta;
									if (allrot[c][ind][dd]>Math.PI)
										allrot[c][ind][dd] -= 2*Math.PI;
									if (allrot[c][ind][dd]<-Math.PI)
										allrot[c][ind][dd] += 2*Math.PI;
								}    
								ind++;
							}
				}

				// ********************************
				// Adding the real rotamer to set.
				// ********************************
				/*    for (int dd=0 ; dd<lib.getChiMax((int) pp[c][2]) ; dd++)
    	allrot[ind][dd] = trueAng[c][dd];
    rmsDevRot[ind] = 0.0;
    probRot[ind] = probRot[pre];
				 */

			} // Of building the library


		// The energy creators that will be used from now on!
		EnergyCreator[] energyCreators1 = {
				new BondCreator(1.0),
				new SimpleHydrogenBond_Dahiyat_LowAccuracy_Creator(1.0),
				new SimpleHPCreator(1.0,1.0,4.25,4.25,false)
		};

		energy = new TotalEnergy(protein, distanceMatrix, energyCreators1, commands);
		hbTerm = (SimpleHydrogenBondEnergy) energy.getEnergyTerm(new SimpleHydrogenBondEnergy());
		solvateTerm = (SimpleHP) energy.getEnergyTerm(new SimpleHP());

		int tmpNum=0,startNum=pp.length-1,endNum=1,inc=-1;
		for (int iter=0 ; iter<TH2 ; iter++) {
			System.out.println("\n*****************************\nDoing Iteration: " + iter + "\n*****************************");	

			double[][] angsAtEndOfIter1 = RotamericTools.getSidechainTorsions(protein);

			if (iter==0)
				w = w1;
			else
				w = w2;

			// Looping on all the residues.
			tmpNum = startNum;
			startNum=endNum;
			endNum = tmpNum;
			if (inc==1)
				inc=-1;
			else
				inc=1;

			for (int s=0 ; s<seder.length ; s++)
				for (int c=0 ; c!=pp.length ; c++)  //for (int c=startNum ; c!=endNum+inc ; c+=inc)  
					if ((pp[c]!=null) && (pp[c][2]==seder[s])) { 
						typ = (int) pp[c][2];
						System.out.println("\n\nDoing residue " + c + " of type " + typ + "\n");	

						inactivateFarFromAtom(protein.residue(c).atoms().findAtomInList("CB", c) , 6.0);


						// ********************************
						// Finding the most similar rotamer
						// ********************************
						pre = RotamericTools.nearestRotFromArray(protein.residue(c),allrot[c]);
						rotlikenative = pre;
						tmprot = allrot[c][pre];
						diffvec[0] = diffvec[1] = diffvec[2] = diffvec[3] = -999;
						for (int dd=0 ; dd<lib.getChiMax((int) pp[c][2]) ; dd++) {
							if (Math.abs(trueAng[c][dd]-tmprot[dd])>Math.PI)
								diffvec[dd] =  180*(2*Math.PI - Math.abs(trueAng[c][dd]-tmprot[dd]))/Math.PI;
							else
								diffvec[dd] =  180*Math.abs(trueAng[c][dd]-tmprot[dd])/Math.PI;
							if (diffvec[dd]>evalLimit*180/Math.PI)
								rotlikenative = -1;
						}   


						//*****************************************************************************    
						// Going over the rotamers
						badPoints = new double[allrot[c].length][];
						finalChis = new double[allrot[c].length][lib.getChiMax((int) pp[c][2])];
						finalRMS = new double[allrot[c].length];
						results = new int[allrot[c].length];
						for (ind=0 ; ind<allrot[c].length ; ind++) {
//							System.out.println("Effetive Rot " + ind);
							for (int dd=0 ; dd<lib.getChiMax((int) pp[c][2]) ; dd++) 
								fourTupple[dd] = allrot[c][ind][dd];

							ResidueBuilder.build(protein.residue(c), 
									protein.residue(c).type, 
									fourTupple);


							/*	energyTerm = energy.getEnergyTerm(new TorsionValEnergy());
	list = ((SimpleEnergyTerm) energyTerm).elementsList();
	for (int dd=0 ; dd<list.size() ; dd++) {
		((TorsionValEnergyElement) list.elementAt(dd)).setTarget(
			fourTupple[((TorsionValEnergyElement) list.elementAt(dd)).getTorCode()]);
//		System.out.println("Set: " + ((TorsionValEnergyElement) list.elementAt(dd)).getTorCode());
	}
	energy.evaluate();
	//minimizer = new LBFGS(energy,0.01,1000,50);
	//System.out.println(minimizer.minimize());

	// Copying the final rotamer to 'finalChis'
    energyTerm = energy.getEnergyTerm(new TorsionValEnergy());
	list = ((SimpleEnergyTerm) energyTerm).elementsList();
	for (int dd=0 ; dd<list.size() ; dd++) 
		finalChis[ind][((TorsionValEnergyElement) list.elementAt(dd)).getTorCode()] =
								((TorsionValEnergyElement) list.elementAt(dd)).torsion().torsion();*/

							energy.evaluate();

							// Copying the final rotamer to 'finalChis'
							for (int dd=0 ; dd<finalChis[ind].length ; dd++) 
								finalChis[ind][dd] = fourTupple[dd];

							// Calculating the RMS
							try {
								finalRMS[ind] = RotamericTools.calcRMS(fullprot.residue(c),protein.residue(c));
							}
							catch (Exception e) {
								finalRMS[ind] = 999;
							}


//							energyTerm = energy.getEnergyTerm(new SolvateEnergy());
//							System.out.println(ind + "\n----------------------");
//							energyTerm.evaluateAtoms();


							// Scores of different energies                           
							badPoints[ind] = brain(c,probRot[c][ind],parametersList); 



							// Checking if the result is true
							diff = Math.abs(finalChis[ind][0]-trueAng[c][0]);
							if (lib.getChiMax((int) pp[c][2])>1) {
								diff1 = Math.abs(finalChis[ind][1]-trueAng[c][1]);
								if (finalChis[ind][1]<0.0)
									diff1alt = Math.abs(finalChis[ind][1]+Math.PI-trueAng[c][1]);
								else
									diff1alt = Math.abs(finalChis[ind][1]-Math.PI-trueAng[c][1]);
							}
							if (Math.abs(diff-2*Math.PI) < diff)
								diff = Math.abs(diff-2*Math.PI);
							if (Math.abs(diff1-2*Math.PI) < diff1)
								diff1 = Math.abs(diff1-2*Math.PI);
							if (Math.abs(diff1alt-2*Math.PI) < diff1alt)
								diff1alt = Math.abs(diff1alt-2*Math.PI);
							if ((pp[c][2]==2) || (pp[c][2]==4) || (pp[c][2]==19))
								diff1 = Math.min(diff1,diff1alt);
							results[ind] = 0;
							if (lib.getChiMax((int) pp[c][2])>1) {
								if ((diff < evalLimit) && (diff1 < evalLimit))
									results[ind] = 1;
								else if (diff < evalLimit)
									results[ind] = 2;
								else if (diff1 < evalLimit)
									results[ind] = 3;
							}
							else {
								if (diff < evalLimit) 
									results[ind] = 1;	
							}

						}  // The loop on the rotamers 
						// END - Going over the rotamers
						//*****************************************************************************    


						// ********************************
						// Making the prediction
						// ********************************
						minscore = 1e10;
						pred = -1;
						for (ind=0 ; ind<allrot[c].length ; ind++) {
							score = w[0][(int) pp[c][2]] * badPoints[ind][0] +
							w[1][(int) pp[c][2]] * badPoints[ind][1] +
							w[2][(int) pp[c][2]] * badPoints[ind][2] +
							w[3][(int) pp[c][2]] * Math.log(badPoints[ind][3]+0.005) +
							w[4][(int) pp[c][2]] * badPoints[ind][4] +
							w[5][(int) pp[c][2]] * badPoints[ind][5];        
							if (score<minscore) {
								pred = ind;
								minscore = score;
							}
						}
						if (pred<0)
							throw new RuntimeException("Serious problems " + pred); 

						// Puting the new prediction
						if ((iter>0) && (w[3][0]>0))
							ResidueBuilder.build(protein.residue(c), 
									protein.residue(c).type, 
									angsAtEndOfIter1[c]);    
						else
							ResidueBuilder.build(protein.residue(c), 
									protein.residue(c).type, 
									allrot[c][pred]);

						// The best rotamer post mortom
						post = 999;
						minrms = 9999;
						for (ind=0 ; ind<allrot[c].length ; ind++) 
							if (finalRMS[ind]<minrms) {
								post = ind;
								minrms = finalRMS[ind];
							}		
						for (int dd=0 ; dd<lib.getChiMax((int) pp[c][2]) ; dd++) {
							preBest[c][dd] = allrot[c][pre][dd];        
							predictBest[c][dd] = allrot[c][pred][dd];        
							postBest[c][dd] = allrot[c][post][dd];
						}
						for (int dd=0 ; dd<lib.getChiMax((int) pp[c][2]) ; dd++) {
							preBestMin[c][dd] = finalChis[pre][dd];        
							predictBestMin[c][dd] = finalChis[pred][dd];        
							postBestMin[c][dd] = finalChis[post][dd];
						}


						// ********************************
						// Writing the raw data
						// ********************************
						System.out.print("\n" + ((iter+1)*111111) + " ");
						System.out.print(typ+ " " + (int) TH1 + " " + c +  "   ");
						System.out.print((rotlikenative+1) + " " + (int)diffvec[0] + " " + (int)diffvec[1] + " " + (int)diffvec[2] + " " + (int)diffvec[3] + "   ");
						System.out.print(fmt1.format(dssp.relACCofRes(c,' ')) + "   ");    
						for (ind=0 ; ind<allrot[c].length ; ind++) 
							System.out.print(
									(ind+1) + " " + fmt3.format(badPoints[ind][0])  + " " + fmt3.format(badPoints[ind][1]) + " " + 
									fmt3.format(badPoints[ind][2])  + " " + fmt2.format(badPoints[ind][3])  + " " +
									fmt.format(badPoints[ind][4])  + " " + fmt.format(badPoints[ind][5])  + "   ");
						for (ind=0 ; ind<allrot[c].length ; ind++) {
							System.out.print(results[ind]  + " ");
							if (ind%3==2)
								System.out.print("    ");
						}
						for (ind=0 ; ind<allrot[c].length ; ind++) {
							System.out.print(fmt1.format(finalRMS[ind])  + " ");
							if (ind%3==2)
								System.out.print("    ");	   
						}
						System.out.println((pre+1) + " " + (pred+1) + " " + (post+1));


					} // Going for the next residue.

			/*    
   // Writing the 6 files to disk
   PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("pre"+(int)TH1+".pdb")));
   for (int kk=0 ; kk<fullprot.residues().size() ; kk++) 
      if (fullprot.residues().residueAt(kk).type > -1) {
      	try {
	    ResidueBuilder.build(fullprot.residues().residueAt(kk), 
	    fullprot.residues().residueAt(kk).type, 
	    preBest[kk]);
            for (int xx=0;xx<fullprot.residues().residueAt(kk).atoms().size();xx++)
                 pw.println(fullprot.residues().residueAt(kk).atoms().atomAt(xx).toString());

	    }
	    catch (Exception e) {
	    	System.out.println("8888 " + e);
	    }
	}
    pw.close();

	pw = new PrintWriter(new BufferedWriter(new FileWriter("preMin"+(int)TH1+".pdb")));
   	for (int kk=0 ; kk<fullprot.residues().size() ; kk++) 
      if (fullprot.residues().residueAt(kk).type > -1) {
      	try {
	    ResidueBuilder.build(fullprot.residues().residueAt(kk), 
	    fullprot.residues().residueAt(kk).type, 
	    preBestMin[kk]);
            for (int xx=0;xx<fullprot.residues().residueAt(kk).atoms().size();xx++)
                 pw.println(fullprot.residues().residueAt(kk).atoms().atomAt(xx).toString());

	    }
	    catch (Exception e) {
	    	System.out.println("8888 " + e);
	    }
	}
    pw.close();

	pw = new PrintWriter(new BufferedWriter(new FileWriter("pred"+(int)TH1+".pdb")));
   	for (int kk=0 ; kk<fullprot.residues().size() ; kk++) 
      if (fullprot.residues().residueAt(kk).type > -1) {
      	try {
	    ResidueBuilder.build(fullprot.residues().residueAt(kk), 
	    fullprot.residues().residueAt(kk).type, 
	    predictBest[kk]);
            for (int xx=0;xx<fullprot.residues().residueAt(kk).atoms().size();xx++)
                 pw.println(fullprot.residues().residueAt(kk).atoms().atomAt(xx).toString());

	    }
	    catch (Exception e) {
	    	System.out.println("8888 " + e);
	    }
	}
    pw.close();

	pw = new PrintWriter(new BufferedWriter(new FileWriter("predMin"+(int)TH1+".pdb")));
   	for (int kk=0 ; kk<fullprot.residues().size() ; kk++) 
      if (fullprot.residues().residueAt(kk).type > -1) {
      	try {
	    ResidueBuilder.build(fullprot.residues().residueAt(kk), 
	    fullprot.residues().residueAt(kk).type, 
	    predictBestMin[kk]);
            for (int xx=0;xx<fullprot.residues().residueAt(kk).atoms().size();xx++)
                 pw.println(fullprot.residues().residueAt(kk).atoms().atomAt(xx).toString());

	    }
	    catch (Exception e) {
	    	System.out.println("8888 " + e);
	    }
	}
    pw.close();

	pw = new PrintWriter(new BufferedWriter(new FileWriter("post"+(int)TH1+".pdb")));
   	for (int kk=0 ; kk<fullprot.residues().size() ; kk++) 
      if (fullprot.residues().residueAt(kk).type > -1) {
      	try {
	    ResidueBuilder.build(fullprot.residues().residueAt(kk), 
	    fullprot.residues().residueAt(kk).type, 
	    postBest[kk]);
            for (int xx=0;xx<fullprot.residues().residueAt(kk).atoms().size();xx++)
                 pw.println(fullprot.residues().residueAt(kk).atoms().atomAt(xx).toString());

	    }
	    catch (Exception e) {
	    	System.out.println("8888 " + e);
	    }
	}
    pw.close();

	pw = new PrintWriter(new BufferedWriter(new FileWriter("postMin"+(int)TH1+".pdb")));
   	for (int kk=0 ; kk<fullprot.residues().size() ; kk++) 
      if (fullprot.residues().residueAt(kk).type > -1) {
      	try {
	    ResidueBuilder.build(fullprot.residues().residueAt(kk), 
	    fullprot.residues().residueAt(kk).type, 
	    postBestMin[kk]);
            for (int xx=0;xx<fullprot.residues().residueAt(kk).atoms().size();xx++)
                 pw.println(fullprot.residues().residueAt(kk).atoms().atomAt(xx).toString());

	    }
	    catch (Exception e) {
	    	System.out.println("8888 " + e);
	    }
	}
    pw.close();
			 */


			// Writing the protein
			pw = new PrintWriter(new BufferedWriter(new FileWriter(writeDir+"pred"+(int)TH1+"_"+iter+".pdb")));
			for (int kk=0 ; kk<protein.residues().size() ; kk++) 
				if ((protein.residues().residueAt(kk).type > -1) && (protein.residues().residueAt(kk).type <20)) {
					for (int xx=0;xx<protein.residues().residueAt(kk).atoms().size();xx++)
						pw.println(protein.residues().residueAt(kk).atoms().atomAt(xx).toString());

				}
			pw.close();

			/*    if (iter>0) {
    	energyTerm = energy.getEnergyTerm(new TorsionValEnergy());
	    list = ((SimpleEnergyTerm) energyTerm).elementsList();
	    for (int dd=0 ; dd<list.size() ; dd++) 
		    ((TorsionValEnergyElement) list.elementAt(dd)).setTarget(
			    predictBest[((TorsionValEnergyElement) list.elementAt(dd)).torsion().getTorsionResNum()]
			    [((TorsionValEnergyElement) list.elementAt(dd)).getTorCode()]);
    	energy.evaluate();
		minimizer = new LBFGS(energy,0.01,10000,50);
		System.out.println(minimizer.minimize());
    	// Writing the minimized protein
		pw = new PrintWriter(new BufferedWriter(new FileWriter(writeDir+"pred"+(int)TH1+"_"+(iter+100)+".pdb")));
   		for (int kk=0 ; kk<protein.residues().size() ; kk++) 
      		if ((protein.residues().residueAt(kk).type > -1) && (protein.residues().residueAt(kk).type <20)) {
            	for (int xx=0;xx<protein.residues().residueAt(kk).atoms().size();xx++)
                	 pw.println(protein.residues().residueAt(kk).atoms().atomAt(xx).toString());
		}
   	 	pw.close();
    }*/





		} // Of the iteration

		System.out.print("\n\n\nSUCCESS IN THE END");        

	}  // of main



	protected static double[] brain(int resnum , double prob , SoftExcludedVolParametersList parametersList) {  
		double tmpenergy;
		SoftExcludedVolParameters parameters;

		double[] result = new double[6];
		for (int c=0 ; c<6 ; c++)
			result[c] = 0;

		result[1] = solvateTerm.hydrophobicEnergy();
		result[2] = solvateTerm.hydrophilicEnergy();
		result[4] = hbTerm.hbEnergy();
		result[5] = -888; 


		result[3] = prob;


		// reseting the relevent atom energies
		AtomList al = protein.residue(resnum).atoms();
		Atom atom,atom1,atom2;
		Iterator atomListIter = al.iterator();
		while ((atom  = (Atom) atomListIter.next()) != null) 
			atom.resetEnergy();
		// Calculating the EV on the relevent atoms  
		Iterator distanceListIter = distanceMatrix.nonBondedList().iterator();
		Distance distance;
		while ((distance  = (Distance) distanceListIter.next()) != null) {
			if (!distance.frozen) {
				atom1 = distance.atom1();
				atom2 = distance.atom2();
				if ((atom1.residueNumber() == resnum) ||
						(atom2.residueNumber() == resnum)) {
					parameters = (SoftExcludedVolParameters) parametersList.parameters(distance);
					if (distance.distance()<parameters.sigma*E) {
						if ((atom1.name().length()==1) || (atom2.name().length()==1) || atom1.name().equals("CA") ||
						atom2.name().equals("CA") || atom1.name().equals("CB") || atom2.name().equals("CB")) {
						tmpenergy = 3*A*(parameters.sigma-distance.distance())*(parameters.sigma-distance.distance());
						}
						else {
						tmpenergy = A*(parameters.sigma-distance.distance())*(parameters.sigma-distance.distance());
						}
//						tmpenergy = A*(parameters.sigma-distance.distance())*(parameters.sigma-distance.distance());
						atom1.addEnergy(tmpenergy/2);
						atom2.addEnergy(tmpenergy/2);
					}
				}
			}
		}
		// Getting the energy values from the sidechain
		atomListIter = al.iterator();
		while ((atom  = (Atom) atomListIter.next()) != null)
			if ((atom.name().length()>1) && 
					(!atom.name().equals("CA")) && 
					(!atom.name().equals("CB")) &&
					(atom.name().charAt(1)!='X')) {
				result[0] += atom.energy();
			}

		return result;
	}

	protected static double[][] readWeight(String fileName) throws Exception {
		double[][] w = new double[6][20];
		String line;
		BufferedReader br = new BufferedReader(new FileReader(fileName));
		for (int c=0 ; c<20 ; c++) {
			line = br.readLine();
			StringTokenizer st = new StringTokenizer(line);
			System.out.println();
			for (int d=0 ; d<6 ; d++) {
				w[d][c] = (new Double(st.nextToken())).doubleValue();
				String tmp = fmt1.format(w[d][c]).toString().trim();
				for (int m=tmp.length() ; m < 5 ; m++)
					tmp += " ";
				System.out.print(tmp + ", ");
			}
		}
		br.close();
		return w;
	}    

	private static void inactivateFarFromAtom(Atom atom , double threshold){
		double threshold2 = threshold*threshold;
		Atom otherAtom;
		AtomList otherResidueAtoms;
		for (int c=0; c<protein.atoms().size() ; c++)
			if (protein.atoms().atomAt(c).active()) 
				protein.atoms().atomAt(c).inActivate();
		for (int c=0; c<protein.atoms().size() ; c++) {
			otherAtom = protein.atoms().atomAt(c);
			if (!otherAtom.active()) 
				if (((atom.x()-otherAtom.x())*(atom.x()-otherAtom.x()) +
				(atom.y()-otherAtom.y())*(atom.y()-otherAtom.y()) +
				(atom.z()-otherAtom.z())*(atom.z()-otherAtom.z())) < threshold2) {
					otherResidueAtoms = protein.residue(otherAtom.residueNumber()).atoms();
					for (int cc=0; cc<otherResidueAtoms.size() ;cc++)
						otherResidueAtoms.atomAt(cc).activate();
				}
		}
	}
	
	/** ================================= init =========================================
	 *
	 *A static function for parsing of the command line arguments and assigning the 
	 *variables commandsFileName, modelFileName and randomNumberSeed with the right inputs. Note that this
	 *static method is using parsing functions such as getOrderedArguments that are defined in MeshiProgram
	 *that MinimizeProtein inherits.
	 **/

	protected static void init(String[] args) {

		/**** NOTE *** the next two lines. Because of a BUG in the Java VM, the 
		 * interfaces "Residues" and "AtomTypes" are not loaded automatically when MinimizeProtein initialize. 
		 * For this purpose these two lines are crucial wherever these two interfaces are implemented. The user might 
		 * rightfully feel that these two lines are "black magic" programming, but happily to our knowledge this is 
		 * the only bizarre phenomenon we are aware of in meshi.
		 **/
		int zvl = ALA; // force the reading of "meshi.parameters.Residues"
		zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"


		String line;
		String errorMessage = ("\n                  ******************\n"+
				"Usage java -Xmx300m NoMinimizePerfect <commands file name> <pdb file name> <index in astral List>\n"+
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));

		commandsFileName = getOrderedArgument(args);
		if (commandsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# commandsFileName = "+commandsFileName);

		modelFileName = getOrderedArgument(args);
		if (modelFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# initial model file name is "+modelFileName);

		dsspFileName = getOrderedArgument(args);
		if (modelFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# dssp file name is "+modelFileName);
		int dirIndex = modelFileName.indexOf('/');
		if (dirIndex==-1)
			writeDir = "";
		else 
			writeDir = modelFileName.substring(0,dirIndex+1);
		System.out.println("# writing prediction to: " + writeDir);

		initRandom(0);

		String tmp = getOrderedArgument(args);
		if (tmp== null) throw new RuntimeException(errorMessage);
		TH1 = (new Integer(tmp)).intValue();
		System.out.println("# The #1 threshhold is:"+TH1);

		tmp = getOrderedArgument(args);
		if (tmp== null) throw new RuntimeException(errorMessage);
		TH2 = (new Integer(tmp)).intValue();
		System.out.println("# The residue type is:"+TH2);

		tmp = getOrderedArgument(args);
		if (tmp== null) throw new RuntimeException(errorMessage);
		A = (new Double(tmp)).doubleValue();
		System.out.println("# A is:"+A);

		tmp = getOrderedArgument(args);
		if (tmp== null) throw new RuntimeException(errorMessage);
		B = (new Double(tmp)).doubleValue();
		System.out.println("# B is:"+B);

		tmp = getOrderedArgument(args);
		if (tmp== null) throw new RuntimeException(errorMessage);
		C = (new Double(tmp)).doubleValue();
		System.out.println("# C is:"+C);

		tmp = getOrderedArgument(args);
		if (tmp== null) throw new RuntimeException(errorMessage);
		D = (new Double(tmp)).doubleValue();
		System.out.println("# D is:"+D);

		tmp = getOrderedArgument(args);
		if (tmp== null) throw new RuntimeException(errorMessage);
		E = (new Double(tmp)).doubleValue();
		System.out.println("# E is:"+E);

		weightsFileName = getOrderedArgument(args);
		if (weightsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# Weights file name is "+weightsFileName);

		weightsFileName1 = getOrderedArgument(args);
		if (weightsFileName1 == null) throw new RuntimeException(errorMessage);
		System.out.println("# Second weights file name is "+weightsFileName1);
	}

} // of class

