package meshi.IMPmy;

public class SA {

	TotalEnergy totalEnergy;
	DomainListWithStructure domainsForReporting;
	
	public SA(TotalEnergy totalEnergy, DomainListWithStructure domainsForReporting) {
		this.totalEnergy = totalEnergy;
		this.domainsForReporting = domainsForReporting;
	}
	
	public void run(int nSteps, double startT, double endT, int reportEvery) {
		int reportNumber = 0;
		for (int step=0 ; step<nSteps ; step++) {
			if (step%reportEvery == 0) {
				domainsForReporting.report(reportNumber);
				reportNumber++;
			}
			
			
			
		}
	}
	
}
