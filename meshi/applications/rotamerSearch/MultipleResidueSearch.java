package meshi.applications.rotamerSearch;

import java.io.IOException;
import java.util.Random;
import java.util.Vector;

import programs.PutHydrogens;

import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.LennardJones.LennardJones;
import meshi.energy.LennardJones.LennardJonesCreator;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBondEnergy;
import meshi.energy.simpleHydrogenBond.SimpleHydrogenBond_Dahiyat_HighAccuracy_Creator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.external.ComplexMESHIconversion;
import meshi.util.file.MeshiWriter;
import meshi.util.rotamericTools.RotamericTools;

public class MultipleResidueSearch   extends MeshiProgram {
	
	public final double maxCummulativeProbInLib = 0.99; // The library will be built up to this cummlative prob, so very rare rotamers are excluded.
	public final double maxLJrisePerRotamer = 15.0; // The maximal rise in LJ energy over base.

	private Protein prot;
	private CommandList commands;
	private Vector<Residue> residues;
	private TotalEnergy energy;
	private DunbrackLib lib;
	private Vector<SingleResidueSearch> singleSearches;
	double[][] pp;
//	LennardJones LJenergy;
//	SimpleHydrogenBondEnergy HBenergy;
	private double baseLJ;
	private double baseHB;
	private double bestScore;
	private int indBestScore;
	private SingleResidueSearch currentlyWorking;
	private double rotamerEnergy;
	
	public MultipleResidueSearch(Protein prot, CommandList commands, Vector<Residue> residues) {
		this.prot = prot;
		this.commands = commands;
		this.residues = residues;
		createEnergyAndRotLib();
		craeteSingleResidueSearch();
	}
	
	private void craeteSingleResidueSearch() {
		int searchOptions = 1;
		singleSearches = new Vector<SingleResidueSearch>();
		for (Residue res : residues) {
			singleSearches.add(new SingleResidueSearch(res, pp[res.number][0], pp[res.number][1], lib));
			searchOptions *= singleSearches.lastElement().howMany();
		}
		System.out.println("Total search options: " + searchOptions);
	}
	
	private void createEnergyAndRotLib() {
		lib  = new DunbrackLib(commands, maxCummulativeProbInLib , 100);
		prot.freeze();
		for (Residue res : residues) {
			res.atoms().defrost();
		}
		EnergyCreator[] energyCreators = {
			new LennardJonesCreator(1.0),
			new SimpleHydrogenBond_Dahiyat_HighAccuracy_Creator(1.0, false)	
		};
		DistanceMatrix distanceMatrix = new DistanceMatrix(prot.atoms(), 5.5, 2.0, 4);  
		pp = RotamericTools.phipsi(prot, distanceMatrix);
		energy = new TotalEnergy(prot, distanceMatrix, energyCreators, commands);
		energy.evaluate();
//		LJenergy = (LennardJones) energy.getEnergyTerm(new LennardJones());
//		HBenergy = (SimpleHydrogenBondEnergy) energy.getEnergyTerm(new SimpleHydrogenBondEnergy());
//		HBenergy.printHBlist();
//		System.exit(0);
		
		baseLJ = energy.energyValues().doubleAt(0);
		baseHB = energy.energyValues().doubleAt(1);
		System.out.print("Initializing exhaustive search for residues: ");
		for (Residue res : residues) {
			System.out.print(res.number + " ");
		}
		System.out.println("\nBase energies: LJ=" + ((int) (baseLJ*100))/100.0 + " HB=" + ((int) (baseHB*100))/100.0 );		
	}
	
	public void search() {
		int counter = 0;
		bestScore = Double.MAX_VALUE;
		currentlyWorking = singleSearches.elementAt(0);
		if (singleSearches.size()==1) {
			for (singleSearches.elementAt(0).intializeSearch() ; singleSearches.elementAt(0).isValid() ; singleSearches.elementAt(0).nextOption()) {
				counter++;
				if (counter%5000 == 0) 
					System.out.println("So far did:" + counter);
				energy.evaluate();
				rotamerEnergy = 
						-(Math.log(0.1+singleSearches.elementAt(0).getProb()));
				if (improvement()) {
					singleSearches.elementAt(0).setBest();
					reportImprovement();
				}				
			}			
			if (bestScore == Double.MAX_VALUE) {
				singleSearches.elementAt(0).restoreInitialCoordinates();
			}
			else {
				singleSearches.elementAt(0).buildBest();
			}
		}
		else if (singleSearches.size()==2) {
			for (singleSearches.elementAt(0).intializeSearch() ; singleSearches.elementAt(0).isValid() ; singleSearches.elementAt(0).nextOption()) {
				for (singleSearches.elementAt(1).intializeSearch() ; singleSearches.elementAt(1).isValid() ; singleSearches.elementAt(1).nextOption()) {
					counter++;
					if (counter%5000 == 0) 
						System.out.println("So far did:" + counter);
					energy.evaluate();
					rotamerEnergy = 
							-(Math.log(0.1+singleSearches.elementAt(0).getProb()) +
							  Math.log(0.1+singleSearches.elementAt(1).getProb()));
					if (improvement()) {
						singleSearches.elementAt(0).setBest();
						singleSearches.elementAt(1).setBest();
						reportImprovement();
					}				
				}
			}			
			if (bestScore == Double.MAX_VALUE) {
				singleSearches.elementAt(0).restoreInitialCoordinates();
				singleSearches.elementAt(1).restoreInitialCoordinates();
			}
			else {
				singleSearches.elementAt(0).buildBest();
				singleSearches.elementAt(1).buildBest();
			}
		}
		else if (singleSearches.size()==3) {
			for (singleSearches.elementAt(0).intializeSearch() ; singleSearches.elementAt(0).isValid() ; singleSearches.elementAt(0).nextOption()) {
				for (singleSearches.elementAt(1).intializeSearch() ; singleSearches.elementAt(1).isValid() ; singleSearches.elementAt(1).nextOption()) {
					for (singleSearches.elementAt(2).intializeSearch() ; singleSearches.elementAt(2).isValid() ; singleSearches.elementAt(2).nextOption()) {
						counter++;
						if (counter%5000 == 0) 
							System.out.println("So far did:" + counter);
						energy.evaluate();
						rotamerEnergy = 
								-(Math.log(0.1+singleSearches.elementAt(0).getProb()) +
										Math.log(0.1+singleSearches.elementAt(1).getProb()) +
										Math.log(0.1+singleSearches.elementAt(2).getProb()));
						if (improvement()) {
							singleSearches.elementAt(0).setBest();
							singleSearches.elementAt(1).setBest();
							singleSearches.elementAt(2).setBest();
							reportImprovement();
						}				
					}
				}			
			}
			if (bestScore == Double.MAX_VALUE) {
				singleSearches.elementAt(0).restoreInitialCoordinates();
				singleSearches.elementAt(1).restoreInitialCoordinates();
				singleSearches.elementAt(2).restoreInitialCoordinates();
			}
			else {
				singleSearches.elementAt(0).buildBest();
				singleSearches.elementAt(1).buildBest();
				singleSearches.elementAt(2).buildBest();
			}
		}
		else {
			throw new RuntimeException("Cannot handle a " + singleSearches.size() + " search at this point.");
		}
	}
	
	/**
	 * If withPermute is false 'r' can be null.
	 */
	public void searchSequential(int Niteration, boolean withPermute , Random r) {
		int[] perm = permutaion(singleSearches.size(),withPermute,r);
		int[] indsBestRot = new int[singleSearches.size()];
		for (int c=0 ; c<indsBestRot.length ; c++) {
			indsBestRot[c] = -1;
		}
		for (int iteration=0 ; iteration<Niteration ; iteration++) {
			System.out.println("Started interation: " + iteration);
			int change = 0;
			for (int c=0 ; c<perm.length ; c++) {
				SingleResidueSearch singleSearch = singleSearches.elementAt(perm[c]);
				currentlyWorking = singleSearch;
				energy.evaluate();
				baseLJ = energy.energyValues().doubleAt(0);
				baseHB = energy.energyValues().doubleAt(1);
				System.out.println("Initializing exhaustive search for residue: " + singleSearch.resName() +  singleSearch.resNumber() + " with possibilities: "  + singleSearch.howMany());
				System.out.println("\nBase energies: LJ=" + ((int) (baseLJ*100))/100.0 + " HB=" + ((int) (baseHB*100))/100.0 );		
				//int counter = 0;
				bestScore = Double.MAX_VALUE;
				for (singleSearch.intializeSearch() ; singleSearch.isValid() ; singleSearch.nextOption()) {
					//counter++;
					//if (counter%5000 == 0) 
					//System.out.println("So far did:" + counter);
					energy.evaluate();
//					System.out.println("Debug: " + )
					rotamerEnergy = 
							-(Math.log(0.1+singleSearch.getProb()));
					if (improvement()) {
						singleSearch.setBest();
						reportImprovement();
					}				
				}
				if (bestScore == Double.MAX_VALUE) {
					singleSearch.restoreInitialCoordinates();
				}
				else {
					singleSearch.buildBest();
				}
				if (indBestScore!=indsBestRot[perm[c]]) {
					indsBestRot[perm[c]] = indBestScore;
					change++;
				}
			}
			System.out.println("Changed " + change + " sidechains.");
		}
	}

	
	/**
	 * This method will evaluate if there is better conformation. Will update bestScore.
	 */
	private boolean improvement() {
		double Elj = energy.energyValues().doubleAt(0);
		double Ehb = energy.energyValues().doubleAt(1);
		if ((baseLJ + maxLJrisePerRotamer *singleSearches.size() ) > Elj) {
			double score = 2.0*Ehb + rotamerEnergy;
			if (score<bestScore) {
				bestScore = score;
				indBestScore = currentlyWorking.pointer();
				return true;
			}
		}
		return false;
	}
	
	private void reportImprovement() {
		System.out.println("\nBest score: " + ((int) (bestScore*100))/100.0 + 
				"   Energies: LJ=" + ((int) (energy.energyValues().doubleAt(0)*100))/100.0 + 
				" HB=" + ((int) (energy.energyValues().doubleAt(1)*100))/100.0 +
				" Rot=" + ((int) (rotamerEnergy*100))/100.0);				
	}
	
	public double bestScore() {
		return bestScore;
	}

	public double indBestScore() {
		return indBestScore;
	}
	
	/**
	 * Permutaion for the sequential treatment
	 */
	public static int[] permutaion(int length, boolean withPermute, Random r) {
		// initialize array and fill it with {0,1,2...}     
		int[] array = new int[length];     
		for(int i = 0; i < array.length; i++)         
			array[i] = i;      
		if (!withPermute) {
			return array;
		}
		else {
			for(int i = 0; i < length; i++){          
				// randomly chosen position in array whose element         
				// will be swapped with the element in position i         
				// note that when i = 0, any position can chosen (0 thru length-1)         
				// when i = 1, only positions 1 thru length -1                     
				// NOTE: r is an insatance of java.util.Random         
				int ran = i + randomNumberGenerator().nextInt(length-i);          
				// perform swap         
				int temp = array[i];         
				array[i] = array[ran];         
				array[ran] = temp;     
			}                            
			return array;
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		initRandom(999);		
//		Protein prot = ComplexMESHIconversion.complex2meshi("C:\\Users\\Nir\\TRiC\\Crystallography\\Intra_ring_interface\\1A6D.pdb");
//		AtomList modelList = new AtomList("C:\\Users\\Nir\\TRiC\\Crystallography\\Intra_ring_interface\\Manual_optimization_on_LJcap_15\\pair_GZ.pdb");
		AtomList modelList = new AtomList("C:\\Users\\Nir\\TRiC\\Crystallography\\Intra_ring_interface\\MESHI_optimization_LJcap_15_ProbBase_10cent\\pair_GZ_iForce.pdb");
		modelList.chainFilter("G").setChain("X");
		modelList.chainFilter("Z").setChain("Y");
		modelList.chainFilter("X").setChain("A");
		modelList.chainFilter("Y").setChain("B");
		Protein prot = ComplexMESHIconversion.complex2meshi(modelList);
		CommandList commands = new CommandList("C:\\Users\\Nir\\MESHI\\commands");
		PutHydrogens.adjustHydrogens(commands, prot);		
		Vector<Residue> residues = new Vector<Residue>();
//		residues.add(prot.residue(387));
		residues.add(prot.residue(1078));
		MultipleResidueSearch search = new MultipleResidueSearch(prot, commands, residues);
		search.search();
//		ComplexMESHIconversion.writeMEHSI2complex(prot, "C:\\Users\\Nir\\TRiC\\Crystallography\\Intra_ring_interface\\1A6D_minimized.pdb");
		Atom.resetNumberOfAtoms();
		AtomList list = ComplexMESHIconversion.MEHSI2complex(prot).duplicate();
		list.chainFilter("A").setChain("#");
		list.chainFilter("B").setChain("^");
		list.chainFilter("#").setChain("G");
		list.chainFilter("^").setChain("Z");
		try {
			list.noOXTFilter().filter(new AtomList.NonHydrogen()).print(new MeshiWriter("C:\\Users\\Nir\\TRiC\\Crystallography\\Intra_ring_interface\\MESHI_optimization_LJcap_15_ProbBase_10cent\\pair_GZ_iter3.pdb"));
//			list.noOXTFilter().filter(new AtomList.NonHydrogen()).print(new MeshiWriter("C:\\Users\\Nir\\TRiC\\Crystallography\\Intra_ring_interface\\Manual_optimization_on_LJcap_15\\pair_BD.pdb"));
		} catch (IOException e) {
			throw new RuntimeException("Could not write file");
		}
	}

}
