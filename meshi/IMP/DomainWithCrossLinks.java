package meshi.IMP;

import java.util.Vector;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.crossLinking.Crosslink;
import meshi.util.crossLinking.CrosslinkVector;

public class DomainWithCrossLinks extends Domain {

	Vector<CrossLinkedLysine> xlLysines;
	Vector<CrossLinkedLysine> boundaries;
	
	public DomainWithCrossLinks(String oneLongString) {
		super(oneLongString);
		xlLysines = null;
		boundaries = null;
	}

		
	public void addCrossLinkData(CrosslinkVector xlVec , AtomList model, boolean toBoundaries) {
		AtomList domainCAatoms = new AtomList();
		for (int c=0 ; c<model.size() ; c++) {
			if (model.atomAt(c).name().equals("CA") && isResNumInDomain(model.atomAt(c).residueNumber())) {
				domainCAatoms.add(model.atomAt(c));
			}
		}		
		if (toBoundaries) {
			boundaries = new Vector<CrossLinkedLysine>();
		}
		else {
			xlLysines = new Vector<CrossLinkedLysine>();			
		}
		for (Crosslink xl : xlVec) {
			if (xl.protName1().equals(proteinName())) {
				if (isResNumInDomain(xl.absPos1())) {
					if (findCrossLinkedLysine(xl.absPos1(),toBoundaries)==null) {
						if (toBoundaries | (domainCAatoms.findAtomInList("CA", xl.absPos1())!=null)) {
							addXLlysine(xl.absPos1(),domainCAatoms,toBoundaries);
						}
					}					
				}				
			}
		}
		for (Crosslink xl : xlVec) {
			if (xl.protName2().equals(proteinName())) {
				if (isResNumInDomain(xl.absPos2())) {
					if (findCrossLinkedLysine(xl.absPos2(),toBoundaries)==null) {
						if (toBoundaries | (domainCAatoms.findAtomInList("CA", xl.absPos2())!=null)) {
							addXLlysine(xl.absPos2(),domainCAatoms,toBoundaries);
						}
					}					
				}				
			}
		}
	}
	
	public CrossLinkedLysine findCrossLinkedLysine(int resNum,boolean toBoundaries) {
		if (toBoundaries) {
			for (CrossLinkedLysine xlK : boundaries) {
				if (xlK.resNum()==resNum) {
					return xlK;
				}
			}
		}
		else {
			for (CrossLinkedLysine xlK : xlLysines) {
				if (xlK.resNum()==resNum) {
					return xlK;
				}
			}			
		}
		return null;
	}

	protected void addXLlysine(int resNum , AtomList model, boolean toBoundaries) {
		double cmx, cmy, cmz; // center of mass x, y and z
		cmx = cmy = cmz = 0.0;
		for (int c=0; c<model.size() ; c++) {
			cmx += model.atomAt(c).x();
			cmy += model.atomAt(c).y();
			cmz += model.atomAt(c).z();
		}
		cmx /= model.size();
		cmy /= model.size();
		cmz /= model.size();
		Atom atom = model.findAtomInList("CA", resNum);
		if (atom==null) {  // Fixing missing atoms in boundaries
			atom = model.findAtomInList("CA", resNum-1); 
		}
		if (atom==null) {  // Fixing missing atoms in boundaries
			atom = model.findAtomInList("CA", resNum+1); 
		}
		if (atom==null) {  // Fixing missing atoms in boundaries
			throw new RuntimeException("Don't know what to do!"); 
		}
		if (toBoundaries) {
			boundaries.add(new CrossLinkedLysine(resNum, (int) (this.x()+(atom.x()-cmx)),
					(int) (this.y()+(atom.y()-cmy)),
					(int) (this.z()+(atom.z()-cmz))));			
		}
		else {
			xlLysines.add(new CrossLinkedLysine(resNum, (int) (this.x()+(atom.x()-cmx)),
					(int) (this.y()+(atom.y()-cmy)),
					(int) (this.z()+(atom.z()-cmz))));
		}
	}
	
	public Vector<CrossLinkedLysine> xlLysines() {
		return xlLysines;
	}

	public Vector<CrossLinkedLysine> boundaries() {
		return boundaries;
	}
}
