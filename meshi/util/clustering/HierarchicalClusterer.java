package meshi.util.clustering;

import java.util.Vector;

import meshi.applications.loopBuilding.AbstractLoopBuilder;

public class HierarchicalClusterer {
	
	private double[][] disMat = null;
	private double[][] largeMatrix = null;	
	private Vector<Cluster> tree = null;
	private int serialTokenIndex = -1;
	private int sizeTokenIndex = -1;
	private int[] sortedByIncreasingSize = null;

	public HierarchicalClusterer (double[][] disMat) {
		this.disMat = disMat;

	}
	
	public void cluster(double maxJoiningDis, int maxClustSize) {
		initialize();
		boolean criteriaNotMet = true;
		for (int step=0; (step<(disMat.length-1)) && criteriaNotMet ; step++) {
			// Finding the minimal distance
			double minDis = Double.MAX_VALUE;
			int minInd1 = -1;
			int minInd2 = -1;
			for (int c1=0 ; c1<tree.size() ; c1++) {
				for (int c2=c1+1 ; c2<tree.size() ; c2++) {
					if ((tree.get(c1).getSonClust()==null) && (tree.get(c2).getSonClust()==null) &&
						((tree.get(c1).getClusterMembers().size()+tree.get(c2).getClusterMembers().size())<=maxClustSize)) {
						if (minDis>largeMatrix[c1][c2]) {
							minDis = largeMatrix[c1][c2];
							minInd1 = c1;
							minInd2 = c2;
						}
					}
				}
			}
			
			// are criteria met?
			if (minDis>maxJoiningDis)
				criteriaNotMet=false;
			else {	// Joining clusters
				//System.out.println("Joining clusters: " + (minInd1) + " " + (minInd2));
				Cluster newClust = new Cluster(largeMatrix);
				newClust.setIndex(tree.size());
				newClust.setClusterR(tree.get(minInd1));
				newClust.setClusterL(tree.get(minInd2));
				tree.get(minInd1).setSonClust(newClust);
				tree.get(minInd2).setSonClust(newClust);
				newClust.setJointDis(minDis);
				newClust.updateClusterMembers();
				// Calculating the distance 
				for (int cc=0 ; cc<tree.size() ; cc++) {
					if (tree.get(cc).getSonClust()==null) {
						double sumDis = 0.0;
						double numDis = 0.0;
						for (int c1=0 ; c1<tree.get(cc).getClusterMembers().size() ; c1++) {
							for (int c2=0 ; c2<newClust.getClusterMembers().size() ; c2++)  {
								sumDis += largeMatrix[tree.get(cc).getClusterMembers().get(c1).intValue()][newClust.getClusterMembers().get(c2).intValue()];
								numDis++;
							}
						}
						largeMatrix[tree.size()][cc] = sumDis/numDis;
						largeMatrix[cc][tree.size()] = sumDis/numDis;
					}
				}
				tree.add(newClust);
			}
		}
	}
	
	public void initializeSerialTokenizer() {
		serialTokenIndex = tree.size();
	}
	
	public Cluster getNextSerialToken() {
		while (serialTokenIndex>0) {
			serialTokenIndex--;
			if (tree.get(serialTokenIndex).getSonClust()==null)
				return tree.get(serialTokenIndex);
		}
		return null;
	}
	
	public void initializeSizeTokenizer() {
		sizeTokenIndex = tree.size();
		// initializing the size array
		double[] tmpSizes = new double[tree.size()];
		for (int c=0 ; c<tree.size() ; c++) {
			tmpSizes[c] = tree.get(c).getClusterMembers().size();
		}
		sortedByIncreasingSize = AbstractLoopBuilder.findTopMinArray(tmpSizes, tmpSizes.length, Double.MAX_VALUE);		
	}
	
	public Cluster getNextSizeToken() {
		while (sizeTokenIndex>0) {
			sizeTokenIndex--;
			if (tree.get(sortedByIncreasingSize[sizeTokenIndex]).getSonClust()==null)
				return tree.get(sortedByIncreasingSize[sizeTokenIndex]);
		}
		return null;
	}

	
	
	private void initialize() {
		tree = new Vector<Cluster>();
		largeMatrix = new double[2*disMat.length-1][2*disMat.length-1];
		for (int i=0 ; i<disMat.length ; i++)
			for (int j=0 ; j<disMat.length ; j++)
				largeMatrix[i][j] = disMat[i][j];
		// make singletons
		for (int i=0 ; i<disMat.length ; i++) {
			tree.add(new Cluster(largeMatrix,i));
		}				
	}
	
	public int getClusterNumber() {
		return tree.size();
	}
	
	public Cluster getClust(int index) {
		return tree.get(index);
	}

}
