package meshi.util.clustering;

public class TestClustering {

	public static void main(String[] args) {
		double[][] p = {{0,2},
				{4,4},
				{1,0},
				{5,5},
				{0,0},
				{0,18}
		};
		
		double[][] mat = new double[p.length][p.length];
		for (int c1=0 ; c1<p.length ; c1++) {
			for (int c2=0 ; c2<p.length ; c2++)  {
				mat[c1][c2] = Math.sqrt((p[c1][0]-p[c2][0])*
						(p[c1][0]-p[c2][0]) +
						(p[c1][1]-p[c2][1])*
						(p[c1][1]-p[c2][1]));
			}
		}
		
		HierarchicalClusterer hc = new HierarchicalClusterer(mat);
		hc.cluster(7.0, 2);
		
	}

}

