package meshi.IMPmy;

import java.io.BufferedWriter;
import java.io.FileWriter;

import meshi.molecularElements.AtomList;
import meshi.util.crossLinking.Crosslink;
import meshi.util.crossLinking.CrosslinkVector;
import meshi.util.file.File2StringArray;

public abstract class DomainListWithStructure extends DomainList {

	private static final long serialVersionUID = 1L;

	public DomainListWithStructure(String fileName) {
		String[] domainText = File2StringArray.f2a(fileName);
		for (int c=0 ; c<domainText.length ; c++) {
			Domain newDomain = new Domain(domainText[c]);
			AtomList model = getAtomicModelFileName(newDomain.proteinName(), newDomain.domainName());
			newDomain.setCenter(model);
			add(newDomain);
		}
	}
	
	public void setXLvecs(CrosslinkVector xlVec, double xlWeight) {
		for (Crosslink xl : xlVec) {
			Domain dom1 = findDomain(xl.protName1(), xl.absPos1());			
			Domain dom2 = findDomain(xl.protName2(), xl.absPos2());
			if (dom1==null) {
				System.out.print("\n\nCould not find domain for residue " + xl.absPos1() + " in protein " + xl.protName1() + "\n\n");
				//throw new RuntimeException("\n\nCould not find domain for residue " + xl.absPos1() + " in protein " + xl.protName1() + "\n\n");
			}
			if (dom2==null) {
				System.out.print("\n\nCould not find domain for residue " + xl.absPos2() + " in protein " + xl.protName2() + "\n\n");
				//throw new RuntimeException("\n\nCould not find domain for residue " + xl.absPos2() + " in protein " + xl.protName2() + "\n\n");
			}
			if ((dom1!=null) & (dom2!=null)) {
				AnchorPosition pos1 = dom1.addPosition(getAtomicModelFileName(dom1.proteinName(),dom1.domainName()),xl.absPos1());
				AnchorPosition pos2 = dom2.addPosition(getAtomicModelFileName(dom2.proteinName(),dom2.domainName()),xl.absPos2());
				disConstList().add(new DistanceConstraint(pos1, pos2, DistanceConstraintType.CROSS_LINK,xlWeight));
			}
		}
	}
	
	public void setBoundaries(double boundaryWeight) {
		boolean handled[] = new boolean[size()];
		for (int c=0 ; c<size() ; c++) {
			handled[c] = false;
		}
		// Main loop on domains
		for (int c=0 ; c<size() ; c++) {
			if (!handled[c]) {
				// getting the scope of the protein;
				String protName = get(c).proteinName();
				int firstRes = get(c).firstRes();
				int lastRes = get(c).lastRes();
				for (int cc=0 ; cc<size() ; cc++) {
					if ((get(cc).proteinName().equals(protName)) && (get(cc).firstRes()<firstRes)) {
						firstRes = get(cc).firstRes();
					}
					if ((get(cc).proteinName().equals(protName)) && (get(cc).lastRes()>lastRes)) {
						lastRes = get(cc).lastRes();
					}
				}
				// Looping on that protein
				for (int res=firstRes+1 ; res<lastRes ; res++) {
					Domain lastDomain = findDomain(protName, res-1);
					Domain thisDomain = findDomain(protName, res);
					if ((lastDomain==null) | (thisDomain==null)) {
						throw new RuntimeException("\n\nCould not find domain for residue " + res + " in protein " + protName + "\n\n");
					}
					if (thisDomain!=lastDomain) {
						if (thisDomain.structured() && lastDomain.structured()) {
						    AnchorPosition pos1 = lastDomain.addPosition(getAtomicModelFileName(lastDomain.proteinName(),lastDomain.domainName()),res-1);
							AnchorPosition pos2 = thisDomain.addPosition(getAtomicModelFileName(thisDomain.proteinName(),thisDomain.domainName()),res);
							disConstList().add(new DistanceConstraint(pos1, pos2, DistanceConstraintType.CONNECTIVITY, boundaryWeight));
						} else {
							// Currently do nothing with unstructured domain but issue a warning
							System.out.print("\n\n******************************************\nNot putting boundary to an unstructred domain.\n******************************************\n\n");
						}
					}	
				}
				// Making sure we do not rehandle that protein
				for (int cc=0 ; cc<size() ; cc++) {
					if (get(cc).proteinName().equals(protName)) {
						handled[cc] = true;
					}
				}
			}
		}
	}
	
	public Domain findDomain(String protName, int resNum) {
		for (Domain domain : this) {
			if (domain.proteinName().equals(protName) && domain.isResNumInDomain(resNum)) {
				return domain;
			}
		}
		return null;
	}
	
	public void putDomainsInRandomCube(double cubeSide) {
		for (Domain dom : this) {
			double newX = cubeSide*(Math.random()-0.5);
			double newY = cubeSide*(Math.random()-0.5);
			double newZ = cubeSide*(Math.random()-0.5);
			dom.moveCenterTo(newX, newY, newZ);
		}
	}

	public void rigidify(double rigidityWeight) {
		for (Domain dom : this) {
			for (int c1=0 ; c1<dom.positions().size() ; c1++) {
				for (int c2=c1+1 ; c2<dom.positions().size() ; c2++) {
					disConstList().add(new DistanceConstraint(dom.positions().get(c1), dom.positions().get(c2), DistanceConstraintType.RIGID_BODY,rigidityWeight));
				}
			}			
			for (int c1=0 ; c1<dom.positions().size() ; c1++) {
				disConstList().add(new DistanceConstraint(dom.positions().get(c1), dom.center(), DistanceConstraintType.RIGID_BODY,rigidityWeight));
			}
		}			
	}

	public void addEV(double evWeight) {
		for (Domain dom : this) {
			dom.updateRtoCenter();
		}
		for (int c1=0 ; c1<size() ; c1++) {
			for (int c2=c1+1 ; c2<size() ; c2++) {
				disConstList().add(new DistanceConstraint(get(c1).center(), get(c2).center(), DistanceConstraintType.EXCLUDED_VOLUME,evWeight));
			}
		}
	}
		
	public void report(int reportNumber) {
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter("report_"+reportNumber+".txt"));

//			boolean handled[] = new boolean[size()];
//			for (int c=0 ; c<size() ; c++) {
//				handled[c] = false;
//			}
//			// Main loop on domains
//			for (int c=0 ; c<size() ; c++) {
//				if (!handled[c]) {
//					// getting the scope of the protein;
//					String protName = get(c).proteinName();
//					int firstRes = get(c).firstRes();
//					int lastRes = get(c).lastRes();
//					for (int cc=0 ; cc<size() ; cc++) {
//						if ((get(cc).proteinName().equals(protName)) && (get(cc).firstRes()<firstRes)) {
//							firstRes = get(cc).firstRes();
//						}
//						if ((get(cc).proteinName().equals(protName)) && (get(cc).lastRes()>lastRes)) {
//							lastRes = get(cc).lastRes();
//						}
//					}
//					// Looping on that protein
//					for (int res=firstRes+1 ; res<lastRes ; res++) {
//						Domain lastDomain = findDomain(protName, res-1);
//						Domain thisDomain = findDomain(protName, res);
//						if ((lastDomain==null) | (thisDomain==null)) {
//							throw new RuntimeException("\n\nCould not find domain for residue " + res + " in protein " + protName + "\n\n");
//						}
//						if (thisDomain!=lastDomain) {
//							if (thisDomain.structured() && lastDomain.structured()) {
//								AnchorPosition pos1 = lastDomain.addPosition(getAtomicModelFileName(lastDomain.proteinName(),lastDomain.domainName()),res-1);
//								AnchorPosition pos2 = thisDomain.addPosition(getAtomicModelFileName(thisDomain.proteinName(),thisDomain.domainName()),res);
//								bw.write(lastDomain.proteinName() + " " + lastDomain.domainName() + " " + pos1.x() + " " + pos1.y() + " " + pos1.z() + "\n");
//								bw.write(thisDomain.proteinName() + " " + thisDomain.domainName() + " " + pos2.x() + " " + pos2.y() + " " + pos2.z() + "\n");
//							} else {
//							}
//						}	
//					}
//					// Making sure we do not rehandle that protein
//					for (int cc=0 ; cc<size() ; cc++) {
//						if (get(cc).proteinName().equals(protName)) {
//							handled[cc] = true;
//						}
//					}
//				}
//			}

//			for (Domain dom : this) {
//				for (AnchorPosition)
//				bw.write(dom.proteinName() + " " + dom.domainName() + " " + dom.center().x() + " " + dom.center().y() + " " + dom.center().z() + "\n");
//			}
			
			for (Domain dom : this) {
				bw.write(dom.proteinName() + " " + dom.domainName() + " " + dom.center().x() + " " + dom.center().y() + " " + dom.center().z() + "\n");
			}
			
			bw.close();
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}  				
	}
	
	abstract AtomList getAtomicModelFileName(String protName, String domainName);
	
}
