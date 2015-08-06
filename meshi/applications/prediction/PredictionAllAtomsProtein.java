package meshi.applications.prediction;

import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.LennardJones.LennardJonesCreator;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.distanceConstrains.TemplateDistanceConstrainsCreator;
import meshi.energy.excludedVol.ExcludedVolCreator;
import meshi.energy.hydrogenBond.HydrogenBondsCreator;
import meshi.energy.hydrogenBondsAngle.HBondsPunishOHNAngleCreator;
import meshi.energy.hydrogenBondsAngle.HbondsPunishHOCAngleCreator;
import meshi.energy.hydrogenBondsPairs.HydrogenBondsPairsCreator;
import meshi.energy.hydrogenBondsPlane.HydrogenBondsPlaneCreator;
import meshi.energy.inflate.InflateCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.propensityTorsion.PropensityTorsionCreator;
import meshi.energy.solvate.SolvateCreatorHBforMinimization;
import meshi.energy.twoTorsions.FlatRamachCreator;
import meshi.energy.twoTorsions.TwoTorsionsCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.ResidueCreator;
import meshi.molecularElements.allAtoms.AllAtomsProtein;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.MinimizerException;
import meshi.parameters.AtomTypes;
import meshi.parameters.MeshiPotential;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.filters.Filter;

public class PredictionAllAtomsProtein  extends AllAtomsProtein implements MeshiPotential, Residues, AtomTypes , KeyWords {
    protected TemplateDistanceConstrainsCreator distanceConstrainsCreator;
    protected Filter originalAtoms;    
    
    public PredictionAllAtomsProtein() {
	super();
    }
    


    public PredictionAllAtomsProtein(AtomList atoms, ResidueCreator creator) {
        super(atoms, creator);
    }

    
    public PredictionLog refine(CommandList commands) throws MinimizerException,LineSearchException {
	System.out.println("---------------------------------- All Atom model refinement -------------------------------------");
	System.out.println("Initial RMS (original residues) = "+
			   PredictionRms.rms(this, commands)+" ("+
			   PredictionRms.rms(this, commands, originalAtoms)+")");
	System.out.println("Initial GDT (original residues) = "+
			   PredictionRms.gdt(this, commands)+" ("+
			   PredictionRms.gdt(this, commands, originalAtoms)+")");
	freeze(new BackboneFilter());
	try {
	     System.out.println("Relaxing with frozen backbone");
	     System.out.println(relax(commands));
	}
	catch (MinimizerException Ex) {
	   System.out.println("Warning: All atom minimization with frozen backbone has failed");
        }
	defrost();
	try {
	     System.out.println("Relaxing with UNFROZEN backbone");
	     System.out.println(relax(commands));
	}
	catch (MinimizerException Ex) {
	   System.out.println("Warning: All atom minimization with unfrozen backbone has failed");
        }
	System.out.println("Initial RMS after short relaxation (original residues) = "+
			   PredictionRms.rms(this, commands)+" ("+
			   PredictionRms.rms(this, commands, originalAtoms)+")");
	System.out.println("Initial GDT after short relaxation (original residues) = "+
			   PredictionRms.gdt(this, commands)+" ("+
			   PredictionRms.gdt(this, commands, originalAtoms)+")");
	AllAtomsMinimizationLoop minimizationLoop = new AllAtomsMinimizationLoop(this, distanceConstrainsCreator, 
							commands);
	try { 
	    minimizationLoop.run();  
	}	
	catch (Exception e) {e.printStackTrace(); throw new RuntimeException("\n"+e+"\n");}
	System.out.println("AfterFirstMinimizationLoop RMS (original residues) = "+
			                   PredictionRms.rms(this, commands)+" ("+
					   PredictionRms.rms(this, commands, originalAtoms)+")");
	System.out.println("AfterFirstMinimizationLoop GDT (original residues) = "+
		            PredictionRms.gdt(this, commands)+" ("+
		            PredictionRms.gdt(this, commands, originalAtoms)+")");
	// run the final minimization with the current conformation as the distance constrain
	//
	distanceConstrainsCreator = new TemplateDistanceConstrainsCreator(originalAtoms);
	distanceConstrainsCreator.setParametersList(this,  commands);

	minimizationLoop = new AllAtomsMinimizationLoop(this, distanceConstrainsCreator, 
							commands);
	minimizationLoop.setTheNumberOfIterationsTo(1);
	try { 
	    minimizationLoop.run();  
	}	

	catch (Exception e) {e.printStackTrace(); throw new RuntimeException("\n"+e+"\n");}
	System.out.println("\nFinal RMS (original residues) = "+
			      PredictionRms.rms(this,commands)+" ("+
			      PredictionRms.rms(this, commands, originalAtoms)+")");

	minimizationLoop.energy().summary();
       	PredictionLog log = new PredictionLog();
	log.comments();
	log.setEnergy(minimizationLoop.energy());
	log.setRms1(PredictionRms.rms(this,commands));
	log.setRms2(PredictionRms.rms(this, commands, originalAtoms));
	log.setGdt1(PredictionRms.gdt(this,commands));
	log.setGdt2(PredictionRms.gdt(this,commands, originalAtoms));
	log.summary1();
	log.summary2();
	log.summary3();
	return log;
    }
   
    //================================= minimization loops =========================================
     private static class AllAtomsMinimizationLoop extends PredictionMinimizationLoop{
	public AllAtomsMinimizationLoop(Protein protein, TemplateDistanceConstrainsCreator distanceConstrainsCreator, 
					CommandList commands) {
	    super(protein, distanceConstrainsCreator, commands, ITERATIONS_ALLATOM);
	}

	public TotalEnergy getEnergy(Protein protein, CommandList commands, 
				     TemplateDistanceConstrainsCreator distanceConstrainsCreator) {
	    HydrogenBondsCreator hbCreator = new HydrogenBondsCreator();
	    EnergyCreator[] energyCreators = {
		new BondCreator(),
		new AngleCreator(),
		new PlaneCreator(),
		new OutOfPlaneCreator(),
		new LennardJonesCreator(),
		hbCreator,
		new HydrogenBondsPairsCreator(hbCreator),
		new HBondsPunishOHNAngleCreator(hbCreator),
		new HbondsPunishHOCAngleCreator(hbCreator),
		new HydrogenBondsPlaneCreator(),
		new SolvateCreatorHBforMinimization(),	                                      
		new FlatRamachCreator(),	                                      
		new TwoTorsionsCreator(),	                                      
		new PropensityTorsionCreator(),
		distanceConstrainsCreator,
		new InflateCreator(),
		new ExcludedVolCreator()
	    };
	    DistanceMatrix distanceMatrix =  new DistanceMatrix(protein.atoms(),5.5,0.5,4); 
	    return new TotalEnergy(protein, distanceMatrix, energyCreators, commands);
	}	
    }



    //================================= filters =========================================
    private static class BackboneFilter implements Filter {
	public boolean accept(Object obj) {
	    if (((Atom) obj).name.equals("CA")) return true;
	    if (((Atom) obj).name.equals("CB")) return true;
	    if (((Atom) obj).name.equals("C")) return true;
	    if (((Atom) obj).name.equals("N")) return true;
	    if (((Atom) obj).name.equals("O")) return true;
	    if (((Atom) obj).name.equals("H")) return true;
	    return false;
	}
    }
    public TemplateDistanceConstrainsCreator distanceConstrainsCreator() {return distanceConstrainsCreator;}
}
