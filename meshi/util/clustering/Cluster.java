package meshi.util.clustering;

import java.util.Arrays;
import java.util.Vector;

public class Cluster {
	
	private double jointDis = -1;
	private Cluster clusterR=null;
	private Cluster clusterL=null;
	private Vector<Integer> clusterMembers = new Vector<Integer>();
	private Cluster sonClust = null;
	private double[][] largeMatrix = null;
	private int index = -1;
	
	public Cluster(double[][] largeMatrix) {
		this.largeMatrix = largeMatrix;
	}

	/**
	 * Make singleton with i
	 */
	public Cluster(double[][] largeMatrix, int i) {
		this.largeMatrix = largeMatrix;
		clusterMembers.add(new Integer(i));
		index = i;
	}

	public double getJointDis() {
		return jointDis;
	}

	public void setJointDis(double jointDis) {
		this.jointDis = jointDis;
	}

	public Cluster getClusterR() {
		return clusterR;
	}

	public void setClusterR(Cluster clusterR) {
		this.clusterR = clusterR;
	}

	public Cluster getClusterL() {
		return clusterL;
	}

	public void setClusterL(Cluster clusterL) {
		this.clusterL = clusterL;
	}

	public Vector<Integer> getClusterMembers() {
		return clusterMembers;
	}

	public void updateClusterMembers() {
		clusterMembers = new Vector<Integer>();
		if (clusterL==null)
			return;
		for (int c=0; c<clusterL.getClusterMembers().size() ; c++) 
			clusterMembers.add(new Integer(clusterL.getClusterMembers().get(c).intValue()));
		for (int c=0; c<clusterR.getClusterMembers().size() ; c++) 
			clusterMembers.add(new Integer(clusterR.getClusterMembers().get(c).intValue()));
	}

	public Cluster getSonClust() {
		return sonClust;
	}

	public void setSonClust(Cluster sonClust) {
		this.sonClust = sonClust;
	}

	public double[][] getLargeMatrix() {
		return largeMatrix;
	}

	public void setLargeMatrix(double[][] largeMatrix) {
		this.largeMatrix = largeMatrix;
	}

	public int getIndex() {
		return index;
	}
	
	public int size() {
		return clusterMembers.size();
	}

	public void setIndex(int index) {
		this.index = index;
	}
	
	public String toString() {
		String str ="";
		str += ("Cluster index: " + index + "\n");
		if (clusterL!=null) 
			str += ("Left parent index: " + clusterL.index + "\n");
		else
			str += ("Left parent index: null\n");
		if (clusterR!=null) 
			str += ("Right parent index: " + clusterR.index + "\n");
		else
			str += ("Right parent index: null\n");
		if (sonClust!=null) 
			str += ("Son index: " + sonClust.index + "\n");
		else
			str += ("Son index: null\n");
		str += ("Number of members: " + clusterMembers.size() + "\n");
		System.out.print(findCenter() + ":   ");
		for (Integer intI : clusterMembers) {
			System.out.print(intI.intValue() + ", ");
		}
		System.out.println();
		return str;
	}

	
	public String toStringWithTrait(double[] trait) {
		String str ="";
		str += ("Cluster index: " + index + "   Size: " + clusterMembers.size() + "   Center: " + findCenter() + "\n");
		str += ("Emin: " + findMinOfTrait(trait) + "   E10: " + findPercentileOfTrait(trait, 0.1) + "   Emean: " + findMeanOfTrait(trait) + "   Ecenter: " + trait[findCenter()] + "\n");
		str += ("-----------------------------------------------------------------------------------------------------\n");
		for (Integer intI : clusterMembers) {
			str += (intI.intValue() + " ");
		}
		str += ("\n\n");
		return str;
	}

	
	/**
	 * NOTICE: Center is found by norm1 distances and not norm2.
	 */
	public int findCenter() {
		double minDis = Double.MAX_VALUE;
		double sumDis = 0.0;
		int minInd = -1;
		for (int cc=0 ; cc<clusterMembers.size() ; cc++) {
			sumDis = 0.0;
			for (int cc1=0 ; cc1<clusterMembers.size() ; cc1++) {
				sumDis += largeMatrix[clusterMembers.get(cc).intValue()][clusterMembers.get(cc1).intValue()];
			}
			if (sumDis<minDis) {
				minDis = sumDis;
				minInd = clusterMembers.get(cc).intValue();
			}
		}
		return minInd;
	}

	public double findMeanOfTrait(double[] trait) {
		double sumEne = 0.0;
		double sumNum = 0.0;
		for (int cc=0 ; cc<clusterMembers.size() ; cc++) {
			sumEne += trait[clusterMembers.get(cc).intValue()];
			sumNum++;
		}
		return sumEne/sumNum;
	}	

	public double findMinOfTrait(double[] trait) {
		double minTrait = Double.MAX_VALUE;
		for (int cc=0 ; cc<clusterMembers.size() ; cc++) {
			if (minTrait > trait[clusterMembers.get(cc).intValue()]) {
				minTrait = trait[clusterMembers.get(cc).intValue()];
			}
		}
		return minTrait;
	}	

	public double findPercentileOfTrait(double[] trait,double prctile) {
		double[] traitToSort = new double[clusterMembers.size()];
		for (int cc=0 ; cc<clusterMembers.size() ; cc++) {
			traitToSort[cc] = trait[clusterMembers.get(cc).intValue()];
		}
		Arrays.sort(traitToSort);
		return traitToSort[(int) (traitToSort.length*prctile)];
	}	
	
}
