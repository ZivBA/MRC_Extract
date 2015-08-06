package meshi.optimizers;
import meshi.energy.TotalEnergy;
import meshi.energy.distanceConstrains.DistanceConstrainsEnergy;
import meshi.energy.distanceConstrains.TemplateDistanceConstrainsCreator;
import meshi.energy.inflate.Inflate;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.Key;
import meshi.util.KeyWords;
                      
public abstract class MinimizationLoop implements KeyWords{
    protected Minimizer minimizer;
    protected TotalEnergy energy;
    protected TemplateDistanceConstrainsCreator distanceConstrainsCreator;
    protected Protein protein;
    protected int nIterations;

    public MinimizationLoop(Protein protein, TemplateDistanceConstrainsCreator distanceConstrainsCreator, CommandList commands , 
			    Key iterationType) {
	nIterations = commands.firstWordFilter(MINIMIZATION_LOOP).secondWord(iterationType).thirdWordInt();
	this.distanceConstrainsCreator = distanceConstrainsCreator;
	this.protein = protein;
	energy = getEnergy(protein, commands, distanceConstrainsCreator);
	minimizer = getMinimizer(energy, commands);
    }

      public MinimizationLoop(Protein protein,
                              TemplateDistanceConstrainsCreator distanceConstrainsCreator,
                              TotalEnergy energy ,
                              CommandList commands,
                              int nIterations) {
	    this.nIterations = nIterations;
            this.distanceConstrainsCreator = distanceConstrainsCreator;
            this.protein = protein;
            this.energy = energy;
            minimizer = getMinimizer(energy, commands);
      }

    
    public abstract TotalEnergy getEnergy(Protein protein, CommandList commands, 
					  TemplateDistanceConstrainsCreator distanceConstrainsCreator);
    public TotalEnergy energy() {return energy;}
					  
    public abstract Minimizer getMinimizer(TotalEnergy energy, CommandList commands);

    public void setTheNumberOfIterationsTo(int nIterations) {
	    	this.nIterations = nIterations;
    }

    public void run() throws MinimizerException, LineSearchException {
	for (int iteration = 1; iteration <= nIterations; iteration++) {
	    System.out.println("\n iteration # "+iteration+"\n");

	    Inflate inflate = (Inflate) energy.getEnergyTerm(new Inflate());
	    DistanceConstrainsEnergy dce = (DistanceConstrainsEnergy) energy.getEnergyTerm(new DistanceConstrainsEnergy());
	    if (inflate == null) throw new RuntimeException("No point in mimizationLoop without inflate");
	    if (iteration != 1) {
		inflate.on();
		if (dce != null) dce.off();
		System.out.println(minimizer.minimize());
	    }
	    inflate.off();
	    dce.on();
	    // here is the simulation
	    try {
		System.out.println(minimizer.minimize());
	    }
	    catch (MinimizerException ex) {
		System.out.println("There is some problem "+ex+"\n"+
				   "I'll try to fix that myself");
		try {
		    inflate.on();
		    if (dce != null) dce.off();
		    System.out.println(minimizer.minimize());
		    inflate.off();
		    dce.on();
		    System.out.println(minimizer.minimize());
		}
		catch (MinimizerException Ex2) {
		    System.out.println("I give up "+Ex2);
		    throw Ex2;
		}
		// Filtering unsatisfied constrains
		if (dce != null) dce.createElementsList(distanceConstrainsCreator.filterParametersList(protein));
	    }
	}
	energy.evaluateAtoms();
   }
    public Protein protein() {return protein;}    
}
	    
	    
	
	
	
    
    
