package programs;

import java.util.Iterator;

import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
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
import meshi.parameters.AtomTypes;
import meshi.parameters.MeshiPotential;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.MeshiProgram;
import meshi.util.rotamericTools.RotamericTools;


/**
 *This class perform concurrent sidechain modeling as describe in Kalisman et al. (2008). There are two modes
 *of operations. One as a stand alone application that prints the new pdb coordinates. Run as follow:
 *java -Xmx300m SCMOD <commands file name> <pdb file name> <number of iterations>
 *
 *The second mode is to give an instance of Protein to the static method "scmod". This mode is useful if one wishes 
 *to do sidechain modeling as part of a larger application. The coordinates of the Protein instance are changed.
 *
 *Important points:
 *-----------------
 *1) The number of iterations can be 0-Inf. If it is 0 all sidechains are put into their most probable back bone 
 *dependent rotamers. 
 *2) In the second mode, residues in the protein that have their CA atoms frozen are not modeled. The stand alone 
 *application models all the residues.
 *3) In the second mode, the protein instance must contain FULL atom residues including all hydrogens.
 *4) In the second mode, you can add a random pertubation to the weights of the rotameric energy. For example, if
 *you put 0.4, then the weight could randomly increase by up to 40%.
 **/


public final class SCMOD extends MeshiProgram implements Residues, 
AtomTypes , KeyWords, MeshiPotential { /**
 * The implemented
 * interfaces defines the 
 * names of atom and residue 
 * types. 
 **/
	private static String commandsFileName = null;
	private static String modelFileName = null; 
	private static int Niter = 1; 

	// The type order in which the side chains are modeled. Starting from TRP (18) and finishing with SER (15)
	private static int[] resTypeModelingOrder = {18,19,4,9,10,12,17,7,6,3,8,14,2,13,11,1,16,15};

	// temporary weight vectors
	private static double[] w1;
	private static double[] w2;
	private static double[] w3;
	private static double[] w4;
	private static double[] w5;
	private static double[] w6;


//	Weights for the first iteration
	// weights for EV
	private static double[] i1w1 = {0.000000, 7.590000, 4.490000, 9.450000, 1.140000, 0.000000, 12.470000, 2.340000, 5.210000, 1.340000, 5.340000, 3.980000, 9.250000, 8.900000, 2.050000, 13.950000, 6.340000, 2.450000, 2.370000, 2.460000};	
	// Weights for hydrophobic-like measure on all the sidechain carbons. (number of solute neighbors - indicating burial)
	private static double[] i1w2 = {0.000000, 0.260000, -0.080000, 0.120000, 0.060000, 0.000000, 0.190000, -0.050000, 0.140000, 0.050000, 0.130000, 0.080000, 0.160000, 0.200000, 0.120000, 0.050000, 0.110000, 0.090000, 0.180000, 0.080000};	
	// Weights for hydrophilic-like measure on all the sidechain polars. (NEGATIVE number of solute neighbors - indicating exposure)
	private static double[] i1w3 = {0.000000, 0.270000, -0.160000, 0.000000, 0.160000, 0.000000, -0.010000, 0.050000, 0.120000, 0.040000, 0.070000, -0.010000, 0.170000, 0.090000, -0.030000, -0.040000, 0.010000, 0.380000, 0.130000, -0.130000};
	// 20 weights for BBDEP rotamer prob
	private static double[] i1w4 = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};	
	// weights for HB tally
	private static double[] i1w5 = {0.000000, 3.680000, 0.520000, 0.290000, 2.200000, 0.000000, 1.530000, 2.900000, 3.280000, 0.370000, 2.790000, 1.480000, 3.130000, 1.560000, 0.860000, 0.390000, 0.600000, 2.560000, 2.970000, 1.040000};	
	// Weights for future use
	private static double[] i1w6 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

//	Weights for iterations 2 till infinity
	// weights for EV
	private static double[] i2w1 = {0.000000, 15.055300, 3.771300, 4.618000, 7.608700, 0.000000, 8.303000, 6.792700, 2.828300, 4.135000, 4.410000, 6.297000, 7.087700, 7.916000, 3.269000, 12.652300, 14.393700, 9.184300, 9.564000, 2.609700};	
	// Weights for hydrophobic like measure on all the carbons. (number of solute neighbors - indicating burial)
	private static double[] i2w2 = {0.000000, 0.293300, 0.111000, 0.012300, 0.152000, 0.000000, 0.475000, 0.132300, 0.055700, 0.052000, 0.160000, 0.091000, 0.101700, 0.124700, 0.087300, 0.054000, 0.135700, 0.193300, 0.290000, 0.078700};	
	// Weights for hydrophilic-like measure on all the sidechain polars. (NEGATIVE number of solute neighbors - indicating exposure)
	private static double[] i2w3 = {0.000000, 0.349700, -0.002000, 0.050300, 0.291000, 0.000000, 0.229000, 0.531300, 0.078300, 0.067000, 0.165000, -0.101000, 0.071300, 0.044700, -0.012000, 0.028000, 0.003700, 0.341000, 0.293000, -0.060700};
	// 20 weights for BBDEP rotamer prob
	private static double[] i2w4 = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};	
	// weights for HB tally
	private static double[] i2w5 = {0.000000, 1.373300, 0.840700, 0.915300, 1.776000, 0.000000, 3.249000, 3.508700, 1.879000, 0.897000, 0.140000, 1.764000, 2.732700, 0.438700, 1.326000, 0.486700, 0.632300, 3.902700, 4.603000, 0.864000};	
	// Weights for future use
	private static double[] i2w6 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


	public static void main(String[] args)  {
		init(args); 
		CommandList commands = new CommandList(commandsFileName);
		Protein protein = new Protein(new AtomList(modelFileName),new ResidueExtendedAtoms(ADD_ATOMS));
		protein.defrost();
		DunbrackLib lib = new DunbrackLib(commands,1.0,100);
		scmod(commands , lib, protein , Niter);
		protein.atoms().print();        
	}


	public static void scmod(CommandList commands , DunbrackLib lib, Protein protein , int maxIter)  {
		scmod(commands , lib, protein , maxIter, 0.0);
	}

	public static void scmod(CommandList commands , DunbrackLib lib, Protein protein , int maxIter, double maxRandomFactor)  {
		scmod(commands , lib, protein , maxIter, maxRandomFactor, null);    	
	}

	public static void scmod(CommandList commands , DunbrackLib lib, Protein protein , int maxIter, double maxRandomFactor, double[][] feedPP)  {
		int pred;
		double score,minscore;
		double[] badPoints;
		double randomFactor = maxRandomFactor*randomNumberGenerator().nextDouble();

		// freezing completely residues that have their CA's frozen.
		for (int c=0; c<protein.residues().size() ; c++) {
			if (protein.residues().residueAt(c).ca() != null) 
				if (protein.residues().residueAt(c).ca().frozen())
					protein.residues().residueAt(c).atoms().freeze();

		}

		SoftExcludedVolParametersList parametersList = 
			new SoftExcludedVolParametersList(commands.firstWord(PARAMETERS_DIRECTORY).secondWord() + "/" + 
					EXCLUDED_VOL_PARAMETERS,1.0);

		DistanceMatrix distanceMatrix = new DistanceMatrix(protein.atoms(), 5.5, 2.0 , 4); 
		double[][] pp;
		if (feedPP != null)
			pp = feedPP;
		else
			pp = RotamericTools.putIntoRot1(protein,distanceMatrix,lib);
		distanceMatrix = new DistanceMatrix(protein.atoms(), 5.5, 2.0 , 4); 

//		********************************
//		Manipulations on the library
//		********************************
		double[][][] allrot = new double[pp.length][][];
		double[][] probRot = new double[pp.length][];
		double[][] rmsDevRot = new double[pp.length][];    
		for (int c=0 ; c!=pp.length ; c++)  
			if ((pp[c]!=null) && (pp[c][2]>0) && (pp[c][2]!=5)) {
				int typ = (int) pp[c][2];
				double rotJump = 15*Math.PI/180.0;
				int eachSide = 0;
				int maxJitterInd1=1,maxJitterInd2=1;
				int ind=0;
				// Jittering the rotamers only for PHE,HIS,TRP and TYR
				if ((typ==4) || (typ==6) || (typ==18) || (typ==19))
					eachSide = 1;
				if ((typ==1) || (typ==15) || (typ==16) || (typ==17)) {
					maxJitterInd1 = 2*eachSide + 1;
					maxJitterInd2 = 1;
				}
				else {
					maxJitterInd1 = 2*eachSide + 1;
					maxJitterInd2 = 2*eachSide + 1;
				}
				int numOfRotamers = lib.getRotamerNum(typ, pp[c][0] , pp[c][1]);
				int effectiveNumOfRotamers = numOfRotamers*maxJitterInd1*maxJitterInd2;
				allrot[c] = new double[effectiveNumOfRotamers][lib.getChiMax(typ)];
				probRot[c] = new double[effectiveNumOfRotamers];
				rmsDevRot[c] = new double[effectiveNumOfRotamers];
				for (int rotInd=0 ; rotInd<numOfRotamers ; rotInd++) {
					double[] tmprot = lib.getRotamer(typ, pp[c][0] , pp[c][1] , rotInd);
					for (int jitterInd1=0 ; jitterInd1<maxJitterInd1 ; jitterInd1++)     
						for (int jitterInd2=0 ; jitterInd2<maxJitterInd2 ; jitterInd2++) {
							probRot[c][ind] = lib.getRotamerProb(typ, pp[c][0] , pp[c][1] , rotInd);
							rmsDevRot[c][ind] = 0.0;
							for (int dd=0 ; dd<lib.getChiMax((int) pp[c][2]) ; dd++) {
								double delta = 0.0;
								if (dd==0)
									delta += (jitterInd1-eachSide)*rotJump;
								if (dd==1)
									delta += (jitterInd2-eachSide)*rotJump;
								rmsDevRot[c][ind] += delta*delta;
								allrot[c][ind][dd] = tmprot[dd] + delta;
								if (allrot[c][ind][dd]>Math.PI)
									allrot[c][ind][dd] -= 2*Math.PI;
								if (allrot[c][ind][dd]<-Math.PI)
									allrot[c][ind][dd] += 2*Math.PI;
							}    
							ind++;
						}
				}
			} // Of building the library



		// The energy creators that will be used from now on!
		EnergyCreator[] energyCreators1 = {
				new SimpleHydrogenBond_Dahiyat_LowAccuracy_Creator(1.0,false),
				new SimpleHPCreator(1.0,1.0,4.25,4.25,false)
		};

		TotalEnergy energy = new TotalEnergy(protein, distanceMatrix, energyCreators1, commands);
		SimpleHydrogenBondEnergy hbTerm = (SimpleHydrogenBondEnergy) energy.getEnergyTerm(new SimpleHydrogenBondEnergy());
		SimpleHP solvateTerm = (SimpleHP) energy.getEnergyTerm(new SimpleHP());

//		The loop of the iterations
		for (int stage=0 ; stage<maxIter ; stage++) {
			System.out.println("\nDoing iteration: " + (stage+1));

			// Setting the weights for this iteration
			if (stage==0) {
				w1 = i1w1;
				w2 = i1w2;    
				w3 = i1w3;    
				w4 = i1w4;    
				w5 = i1w5;    
				w6 = i1w6;
			}
			else {
				w1 = i2w1;
				w2 = i2w2;    
				w3 = i2w3;    
				w4 = i2w4;    
				w5 = i2w5;    
				w6 = i2w6;
			}    	

			// Applying the random modification to the rotameric energies
			for (int c=0 ; c<w4.length ; c++) {
				w4[c] += w4[c]*randomFactor;
			}


			// Looping on the type order
			for (int cc=0 ; cc<resTypeModelingOrder.length ; cc++) {
				// Looping on the chain from N to C termini
				for (int c=0 ; c<pp.length ; c++) {
					// Valid residue check:	It exist, It is of the type we currently model, It is not frozen
					if ((pp[c]!=null) && (pp[c][2]==resTypeModelingOrder[cc]) && !protein.residue(c).ca().frozen()) {
						int typ = (int) pp[c][2];
						minscore = 1e10;
						pred = -1;
						inactivateFarFromAtom(protein, protein.residue(c).atoms().findAtomInList("CB", c) , 6.0);

						//*****************************************************************************    
						// Going over the rotamers
						for (int ind=0 ; ind<allrot[c].length ; ind++) {
							// puting the rotamer in the protein
							ResidueBuilder.build(protein.residue(c), 
									protein.residue(c).type, 
									allrot[c][ind]);
							// Calculating the energies for this rotamer
							energy.evaluate();
							badPoints = calcVariousEnergies(c,probRot[c][ind],parametersList, solvateTerm, hbTerm, protein, distanceMatrix);
							score = w1[typ]*badPoints[0] +
							w2[typ]*badPoints[1] + 
							w3[typ]*badPoints[2] +
							w4[typ]*badPoints[3] +
							w5[typ]*badPoints[4] +
							w6[typ]*badPoints[5];  
							if (score<minscore) {
								pred = ind;
								minscore = score;
							}
						}  // The loop on the rotamers 
						// END - Going over the rotamers
						//*****************************************************************************    

						// Putting the prediction in the protein
						if (pred<0)
							throw new RuntimeException("Serious problems " + pred); 
						ResidueBuilder.build(protein.residue(c),typ,allrot[c][pred]);

					} // This is a valid residue check
				} // Going for the next residue.
			} // Going for the next residue type.
		} // Going for the next iteration.

	}  // scmod


	protected static double[] calcVariousEnergies(int resnum , double prob , SoftExcludedVolParametersList parametersList, 
			SimpleHP solvateTerm, SimpleHydrogenBondEnergy hbTerm, Protein protein, DistanceMatrix distanceMatrix) {  
		double tmpenergy;
		SoftExcludedVolParameters parameters;

		double[] result = new double[6];
		for (int c=0 ; c<6 ; c++)
			result[c] = 0;

		result[1] = solvateTerm.hydrophobicEnergy();
		result[2] = solvateTerm.hydrophilicEnergy();
		result[4] = hbTerm.hbEnergy();


		result[3] = Math.log(prob + 0.005);


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
					if (distance.distance()<parameters.sigma) {
						if ((atom1.name().length()==1) || (atom2.name().length()==1) || atom1.name().equals("CA") ||
						atom2.name().equals("CA") || atom1.name().equals("CB") || atom2.name().equals("CB")) {
						tmpenergy = 3*(parameters.sigma-distance.distance())*(parameters.sigma-distance.distance());
						}
						else {
						tmpenergy = (parameters.sigma-distance.distance())*(parameters.sigma-distance.distance());
						}
//						tmpenergy = (parameters.sigma-distance.distance())*(parameters.sigma-distance.distance());
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


	/**
	 * For faster run times, I inactivate atoms that are far from the modeled sidechain.
	 */
	private static void inactivateFarFromAtom(Protein protein, Atom atom , double threshold){
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


		String errorMessage = ("\n                  ******************\n"+
				"Usage java -Xmx300m SCMOD <commands file name> <pdb file name> <number of iterations>\n"+
		"                    ******************\n");

		if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));

		commandsFileName = getOrderedArgument(args);
		if (commandsFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# commandsFileName = "+commandsFileName);

		modelFileName = getOrderedArgument(args);
		if (modelFileName == null) throw new RuntimeException(errorMessage);
		System.out.println("# initial model file name is "+modelFileName);

		String iterString = getOrderedArgument(args);
		if (iterString== null) throw new RuntimeException(errorMessage);
		Niter = (new Integer(iterString)).intValue();
		System.out.println("# Number of iterations "+ Niter);


		initRandom(0);

	}

} // of class

