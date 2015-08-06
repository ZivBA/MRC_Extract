package meshi.IMPmy;

import java.util.Vector;

import meshi.util.crossLinking.CrosslinkVector;

public abstract class DomainList  extends Vector<Domain>  {
	
	private static final long serialVersionUID = 1L;

	private DistanceConstraintList disConstList;
	
	public DomainList() { 
		disConstList = new DistanceConstraintList();
	}
	
	public abstract void setXLvecs(CrosslinkVector xlVec, double xlWeight);
	
	public abstract void setBoundaries(double boundaryWeight);
	
	public abstract void rigidify(double rigidityWeight); 

	public abstract void addEV(double evWeight);
		
	public abstract void report(int reportNumber); 
	
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

	public DistanceConstraintList disConstList() {
		return disConstList;
	}
	
	public void resetConstraintList() {
		disConstList = new DistanceConstraintList();
	}

}
