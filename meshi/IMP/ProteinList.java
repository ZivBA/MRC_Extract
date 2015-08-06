package meshi.IMP;

import java.util.Vector;

import meshi.applications.TriC.TricAlignment;
import meshi.molecularElements.AtomList;
import meshi.util.crossLinking.Crosslink;
import meshi.util.crossLinking.CrosslinkVector;
import meshi.util.crossLinking.CrosslinkVectorTRiC;
import meshi.util.file.File2StringArray;

public class ProteinList extends Vector<Protein> {

	public ProteinList() {	}

	public ProteinList(String fileName) {
		String[] domainText = File2StringArray.f2a(fileName);
		for (int c=0 ; c<domainText.length ; c++) {
			DomainWithCrossLinks newDomain = new DomainWithCrossLinks(domainText[c]);
			Protein prot = findProteinByName(newDomain.proteinName());
			if (prot==null) {
				prot = new Protein(newDomain.proteinName());
				add(prot);
			}
			prot.addDomain(newDomain);			
		}		
	}
	
	public Protein findProteinByName(String name) {
		for (Protein protein : this) {
			if (protein.proteinName().equals(name)) {
				return protein;
			}
		}
		return null;
	}
	
	public String toString() {
		String str = "";
		for (Protein prot : this) {
			str += prot.toString();
		}
		return str;
	}
	
	public DomainList allDomains() {
		DomainList allDomains = new DomainList();
		for (Protein protein : this) {
			for (DomainWithCrossLinks domain : protein.domainList()) {
				allDomains.add(domain);	
			}
		}
		return allDomains;
	}
	
	public void putDomainsRandomlyInSphere(double initialRadius) {
		DomainList allDomains = allDomains();
		for (Domain domain : allDomains) {
			int radius = (int) (initialRadius*Math.random());
			int x = (int) (1000*(Math.random()-0.5));
			int y = (int) (1000*(Math.random()-0.5));
			int z = (int) (1000*(Math.random()-0.5));
			int currentRadius = (int) Math.sqrt(x*x+y*y+z*z);
			x = (x*radius/currentRadius);
			y = (y*radius/currentRadius);
			z = (z*radius/currentRadius);			
			domain.setXYZ(x, y, z);			
		}		
	}
	
	
	public CrosslinkVector getBoundariesAsCrossLinks() {
		CrosslinkVectorTRiC vec = new CrosslinkVectorTRiC(); // No one should know it's TRiC
		for (Protein prot : this) {
			DomainList domList = prot.domainList();
			int firstRes = domList.elementAt(0).firstRes();
			int lastRes = domList.elementAt(0).lastRes();
			for (Domain domain : domList) {
				if (domain.lastRes()>lastRes) {
					lastRes = domain.lastRes();
				}
			}
			for (int resNum=firstRes ; resNum<lastRes ; resNum++) {
				Domain thisDom = prot.findDomain(resNum);
				Domain nextDom = prot.findDomain(resNum+1);
				if ((thisDom!=null) && (nextDom!=null) && (thisDom!=nextDom)) {
					Crosslink xl = new Crosslink(-1,
							"",
							"", 
							-1, 
							-1, 
							prot.proteinName(),
							prot.proteinName(),
							"",
							resNum, 
							resNum+1);
					vec.add(xl);
				}
			}
		}
		return vec;
	}


	public void addCrossLinkAndBoundaryData(CrosslinkVector xlVec) {
		DomainList allDomains = allDomains();
		CrosslinkVector boundaryList = getBoundariesAsCrossLinks();
		for (DomainWithCrossLinks domain : allDomains) {
			domain.addCrossLinkData(xlVec, new AtomList("HM_"+domain.proteinName()+".pdb"),false);
			domain.addCrossLinkData(boundaryList, new AtomList("HM_"+domain.proteinName()+".pdb"),true);
		}		
	}

		
	
	public static void main(String[] args) {
		ProteinList protList = new ProteinList("domainsBovine.txt");
		protList.putDomainsRandomlyInSphere(300.0);
		CrosslinkVector xlVec = new CrosslinkVectorTRiC("Abersold_Bovine_intra_inter.txt",0).filterOutIntraUnit();
		CrosslinkVector xlVec1 = new CrosslinkVectorTRiC("NK_combined_moderate.txt",1).filterOutIntraUnit();
		xlVec.addVec(xlVec1);
		protList.addCrossLinkAndBoundaryData(xlVec);
		RepresentationIMP rep = new RepresentationIMP(protList.allDomains());
		rep.writeToDisk("TRiC_1_ring");
		RestraintIMP rest = new RestraintIMP(protList.allDomains());
		rest.writeToDisk("TRiC_1_ring", xlVec,protList.getBoundariesAsCrossLinks());
		DisplayIMP disp = new DisplayIMP(protList.allDomains());
		disp.writeToDisk("TRiC_1_ring");
		
		
		

		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
//		String[] unitNames = {"A",
//				"B",
//				"G",
//				"D",
//				"E",
//				"H",
//				"Q",
//				"Z"};
//		
//		int[][][] domainsOn1Q3R = {{{15,149} , {404,526}},
//				{{150,217} , {369,403}},
//				{{218,247} , {283,368}}};
//		
//		TricAlignment alignment = new TricAlignment();
//		
//		for (int unit=0 ; unit < unitNames.length ; unit++) {
//			for (int dom=0 ; dom<3 ; dom++) {
//				System.out.print(unitNames[unit]+" "+dom);
//				for (int part=0 ; part<2 ; part++) {
//					System.out.print(" " + alignment.getNewResNum('K', domainsOn1Q3R[dom][part][0], unitNames[unit].charAt(0)));
//					System.out.print(" " + alignment.getNewResNum('K', domainsOn1Q3R[dom][part][1], unitNames[unit].charAt(0)));					
//				}
//				System.out.println();
//			}
//		}
		
		
	}
	
	
}
