package meshi.IMPmy;

public class SteepestDescent {

	TotalEnergy totalEnergy;
	DomainList domainsForReporting;
	
	public SteepestDescent(TotalEnergy totalEnergy, DomainList domainsForReporting) {
		this.totalEnergy = totalEnergy;
		this.domainsForReporting = domainsForReporting;
	}

	// Returns 0 if converged. 1 if not. 
	public int run(int maxSteps, double gradTol, int reportEvery) {
		int stepCounter = 0;
		int reportNumber = 0;
		double energyInLastPosition = totalEnergy.evaluate(true);
		double energyInNewPosition = 0.0;
		double gradMagInLastPosition = totalEnergy.getGradMagnitude();
		double gradMult =  Math.min(1.0 , 2.0/gradMagInLastPosition);
		while ((stepCounter<maxSteps) && (gradMagInLastPosition>gradTol)) {
			if (stepCounter%reportEvery == 0) {
				domainsForReporting.report(reportNumber);
				System.out.println("Minimization step: " + stepCounter + "   Energy: " + energyInLastPosition + "   Gradient: " + gradMagInLastPosition + "   GradMult: " + gradMult);
				reportNumber++;
			}
			totalEnergy.snapshot();
			energyInNewPosition = Double.MAX_VALUE;
			boolean successfulStep = false;
			while (!successfulStep) {
				totalEnergy.moveGrad(gradMult);
				energyInNewPosition = totalEnergy.evaluate(false);
				if (energyInNewPosition>energyInLastPosition) {
					totalEnergy.restoreSnapshot();
					gradMult *= 0.5;
				}
				else {
					successfulStep = true;
				}
			}
			energyInLastPosition = totalEnergy.evaluate(true);
			gradMagInLastPosition = totalEnergy.getGradMagnitude();
			gradMult *= 1.5;
			stepCounter++;
		}
		if (stepCounter==maxSteps) {
			return 1;
		}
		else {
			return 0;
		}
	}
}
