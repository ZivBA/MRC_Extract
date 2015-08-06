package meshi.util.crossLinking;

public class Compare2XLsets {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		CrosslinkVectorPIC vec1 = new CrosslinkVectorPIC(args[0], 1);
		CrosslinkVectorPIC vec2 = new CrosslinkVectorPIC(args[1], 1);
		System.out.println(vec1.compareToOtherXLvec(vec2));
	}

}
