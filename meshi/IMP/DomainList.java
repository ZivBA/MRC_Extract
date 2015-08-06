package meshi.IMP;

import java.util.Vector;

public class DomainList extends Vector<DomainWithCrossLinks> {
	
	public DomainList() {}
	
	public Domain findDomainByName(String name) {
		for (Domain domain : this) {
			if (domain.domainName().equals(name)) {
				return domain;
			}
		}
		return null;
	}
	
	public Domain findDomainByResidue(int resNum) {
		for (Domain domain : this) {
			for (int c=0 ; c<domain.domainResidues().length ; c++) {
				if ((domain.domainResidues()[c][0]<=resNum) &&
						(domain.domainResidues()[c][1]>=resNum)) {
					return domain;
				}
			}
		}
		return null;
	}

	public Domain findDomain(String protName, int resNum) {
		for (Domain domain : this) {
			if (domain.proteinName().equals(protName)) {
				for (int c=0 ; c<domain.domainResidues().length ; c++) {
					if ((domain.domainResidues()[c][0]<=resNum) &&
							(domain.domainResidues()[c][1]>=resNum)) {
						return domain;
					}
				}
			}
		}
		return null;
	}
		
	public String toString() {
		String str = "";
		for (Domain domain : this) {
			str += domain.toString();
		}
		return str;
	}


}
