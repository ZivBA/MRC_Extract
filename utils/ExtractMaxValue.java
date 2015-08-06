package utils;
/**
 * Created by Ziv_BA on 28/07/2015.
 */
public class ExtractMaxValue {
	public static float[] getMaxValue(MRC_Map_New mrcMap){

		float[][][] grid = mrcMap.getGrid();
		// get the dimentions of the map file (x,y,z max values)
		int cordX = grid.length;
		int cordY = grid[1].length;
		int cordZ = grid[1][1].length;
		// get the grid file
		//TODO possibly improve code by removing cordXYZ and using the grid.
		float maxIntensitiValue = Float.MIN_VALUE;
		float[] maxCoords = new float[0];
		for (int x = 0; x < cordX; x++) {
			for (int y = 0; y < cordY; y++) {
				for (int z = 0; z < cordZ; z++) {
					float tempValue = grid[x][y][z];
					if (tempValue > maxIntensitiValue){
						maxIntensitiValue = tempValue;
						maxCoords = mrcMap.getActualCoords(x,y,z);
					}
				}
			}
		}
		return new float[] {maxIntensitiValue, maxCoords[0],maxCoords[1],maxCoords[2]};
	}
}
