package meshi.optimizers;
import java.util.Random;

import meshi.energy.TotalEnergy;
import meshi.energy.inflate.Inflate;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.MeshiProgram;
                      
public class MCM implements KeyWords{
    protected Minimizer minimizer;
    protected TotalEnergy energy;
    protected Protein protein;
    protected int nIterations;
    protected double initialTemperature;
    protected double finalTemperature;
    protected double deltaT;
    public MCM(Protein protein, TotalEnergy energy, CommandList commands) {
	nIterations = commands.firstWordFilter(MCM).secondWord(N_MCM_STEPS).thirdWordInt();
	initialTemperature = commands.firstWordFilter(MCM).secondWord(INITIAL_TEMPERATURE).thirdWordDouble();
	finalTemperature = commands.firstWordFilter(MCM).secondWord(FINAL_TEMPERATURE).thirdWordDouble();
	deltaT = (finalTemperature - initialTemperature)/(nIterations+1);
	this.energy = energy;
	minimizer = new LBFGS(energy, commands);
    }


    public void run() throws MinimizerException, LineSearchException {
	double oldEnergy = 1000000;
	double[][] oldCoordinates;
	double currentEnergy;
	double dE;
	Random randomNumberGenerator = MeshiProgram.randomNumberGenerator();

	for (int iteration = 1; iteration <= nIterations; iteration++) {
	    double temperature = initialTemperature +iteration*deltaT;
	    System.out.println("\n iteration # "+iteration+"; temperature "+temperature+" current energy = "+oldEnergy+"\n");
	    oldCoordinates = getOldCoordinates(energy);	    
	    Inflate inflate = (Inflate) energy.getEnergyTerm(new Inflate());
	    if (inflate == null) throw new RuntimeException("No point in MCM without inflate");
	    if (iteration != 1) {
		inflate.on();
		System.out.println(minimizer.minimize());
	    }
	    inflate.off();
	    System.out.println(minimizer.minimize());
	    currentEnergy = energy.evaluate();
	    dE = currentEnergy - oldEnergy;
	    if (dE > 0) {
		double rnd = randomNumberGenerator.nextDouble();
		if (rnd > Math.exp(-dE/temperature)) {//That is Metropolis criterion failed
		    setCoordinates(energy,oldCoordinates);
		    System.out.println("This step failed: dE = "+dE+
				       "  dE/temperature = "+ dE/temperature+
				       "  Math.exp(-dE/temperature) = "+Math.exp(-dE/temperature)+
				       "  rnd = "+rnd); 
		}
		else  {
		    System.out.println("This step succeeded: dE = "+dE+
				       "  dE/temperature = "+ dE/temperature+
				       "  Math.exp(-dE/temperature) = "+Math.exp(-dE/temperature)+
				       "  rnd = "+rnd); 
		    oldEnergy = currentEnergy;
		}
	    }
	    else {
		System.out.println("This step succeeded: dE = "+dE);
		oldEnergy = currentEnergy;
	    }
	}
    }
    
    private static double[][] getOldCoordinates(TotalEnergy energy) {
	double[][] coordinates = energy.coordinates();
	int length = coordinates.length;
	double[][] out = new double[length][2];

	for (int i = 0; i < length; i++) {
	    out[i][0] = coordinates[i][0];
	    out[i][1] = coordinates[i][1];
	}
	
	return out;
    }

    private static void setCoordinates(TotalEnergy energy, double[][] toSet) {
	double[][] coordinates = energy.coordinates();
	int length = coordinates.length;
	if (length != toSet.length) throw new RuntimeException("Weird parameters to MCM.setCoordinates");

	for (int i = 0; i < length; i++) {
	    coordinates[i][0] = toSet[i][0];
	    coordinates[i][1] = toSet[i][1];
	}
    }
}    
	    
	    
	
	
	
    
