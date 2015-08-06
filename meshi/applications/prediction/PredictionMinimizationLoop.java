package meshi.applications.prediction;
import meshi.energy.TotalEnergy;
import meshi.energy.distanceConstrains.TemplateDistanceConstrainsCreator;
import meshi.molecularElements.Protein;
import meshi.optimizers.LBFGS;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.MinimizationLoop;
import meshi.optimizers.Minimizer;
import meshi.optimizers.MinimizerException;
import meshi.util.CommandList;
import meshi.util.Key;

public abstract class PredictionMinimizationLoop extends MinimizationLoop{
    public PredictionMinimizationLoop(Protein protein, 
			    TemplateDistanceConstrainsCreator distanceConstrainsCreator, 
			    CommandList commands , Key iterationType) {
	super(protein, distanceConstrainsCreator, commands , iterationType);
    }

    public Minimizer getMinimizer(TotalEnergy energy, CommandList commands){
	return new LBFGS(energy, commands);
    }
    public void run() throws MinimizerException, LineSearchException {
	super.run();
    }
}
    
	    
