package meshi.IMP;

public class Protein {

	private String proteinName;
	private DomainList domainList;
	
	public Protein(String proteinName, DomainList domainList) {
		this.proteinName =  proteinName;
		this.domainList = domainList;
	}
	
	public Protein(String protName) {
		this(protName, new DomainList());
	}
	
	public void addDomain(DomainWithCrossLinks newDomain) {
		domainList.add(newDomain);
	}

	public String proteinName() {
		return proteinName;
	}
	
	public String toString() {
		return domainList.toString();
	}

	public DomainList domainList() {
		return domainList;
	}
	
	public Domain findDomain(int resNum) {
		for (Domain domain : domainList) {
			if (domain.isResNumInDomain(resNum)) {
				return domain;
			}
		}
		return null;
	}

}
