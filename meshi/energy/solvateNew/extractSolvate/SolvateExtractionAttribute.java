package meshi.energy.solvateNew.extractSolvate;

import meshi.util.MeshiAttribute;

public class SolvateExtractionAttribute implements MeshiAttribute {
	
	private double sumOfNeighbors = 0.0;
	private double sumOfHB = 0.0;
	private double sumOfSB = 0.0;
	
	public SolvateExtractionAttribute() {}

	public int key() {
		return SOLVATE_EXTRACTION_ATTRIBUTE;
	}

	public double getSumOfNeighbors() {
		return sumOfNeighbors;
	}

	public void addToSumOfNeighbors(double add) {
		sumOfNeighbors += add;
	}

	public double getSumOfHB() {
		return sumOfHB;
	}

	public void addToSumOfHB(double add) {
		sumOfHB += add;
	}

	public double getSumOfSB() {
		return sumOfSB;
	}

	public void addToSumOfSB(double add) {
		sumOfSB += add;
	}
	

}
