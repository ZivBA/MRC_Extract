package utils.ScoreUtilities;

import utils.UtilExceptions.CoordOutOfRangeException;

import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
public class MRC_Map_New {


	protected int Nx, Ny, Nz;
	protected int dataType;
	protected int startX, startY, startZ;
	protected int intervalX, intervalY, intervalZ;
	protected float dimX, dimY, dimZ;
	protected float angX, angY, angZ;
	protected float Dmin, Dmax, Dmean;
	protected float Ox, Oy, Oz;
	protected float[][][] grid;
	protected float[] xAxis, yAxis, zAxis;

	public MRC_Map_New(String MRC_filename) {
		System.out.print("Reading map: " + MRC_filename + "\n");
		readMRC(MRC_filename);
		if (MRC_filename.endsWith("piciis1.mrc") || MRC_filename.endsWith("pic_iis1_flip_ok_pixel.mrc")
				|| MRC_filename.endsWith("pic_del_tri1_flip_ok_pixel.mrc") || MRC_filename.endsWith("E_and_H.mrc")
				|| MRC_filename.endsWith("PolII_map_portion.mrc")) {
			System.out.println("\n\nWarning!!! Warning!!! Warning!!!    Changing the origin of the map!\n");
			Ox = (float) -185.0;
			Oy = (float) -185.0;
			Oz = (float) -185.0;
		}
		if (MRC_filename.endsWith("My_delineated_EH_factors.mrc")) {
			System.out.println("\n\n Warning!!! Warning!!! Warning!!!  The Origins should be [0,0,0], but it's only" +
					" close to that [" + Ox + "," + Oy + "," + Oz + "]. Changing to [-0.01,-0.01,-0.01]\n");
			Ox = (float) -0.01;
			Oy = (float) -0.01;
			Oz = (float) -0.01;
		}
		xAxis = new float[Nx];
		float v = Ox;
		for (int ix = 0; ix < Nx; ix++) {
			xAxis[ix] = v;
			v += (dimX / Nx);
		}
		yAxis = new float[Ny];
		v = Oy;
		for (int iy = 0; iy < Ny; iy++) {
			yAxis[iy] = v;
			v += (dimY / Ny);
		}
		zAxis = new float[Nz];
		v = Oz;
		for (int iz = 0; iz < Nz; iz++) {
			zAxis[iz] = v;
			v += (dimZ / Nz);
		}
		System.out.println("The grid size is: [" + Nx + " " + Ny + " " + Nz + "]");
		System.out.println("Data type is: " + dataType);
		System.out.println("The start indices are: [" + startX + " " + startY + " " + startZ + "]");
		System.out.println("The intervals in each side are: [" + intervalX + " " + intervalY + " " + intervalZ + "]");
		System.out.println("The cell dimensions are: [" + dimX + " " + dimY + " " + dimZ + "]");
		System.out.println("The cell angles are (currently can only handle boxes): [" + angX + " " + angY + " " + angZ + "]");
		System.out.println("Minimum intensity value: " + Dmin + "   Maximum intensity value: " + Dmax + "   Mean " +
				"intensity value: " + Dmean);
		System.out.println("The cell origin is: [" + Ox + " " + Oy + " " + Oz + "]");
		System.out.println("Map read successfully.");

	}


	public float[][][] getGrid() {
		return grid;
	}

	public double val(double x, double y, double z) {
		if ((x > xAxis[xAxis.length - 1]) |
				(y > yAxis[yAxis.length - 1]) | // was "y > xAxis[..]", fixed.
				(z > zAxis[zAxis.length - 1]) |
				(x < xAxis[0]) |
				(y < yAxis[0]) |
				(z < zAxis[0])) {
//			return -0.1*Math.sqrt(x*x+y*y+z*z);
			throw new CoordOutOfRangeException("\nIndex out of range: [" + x + "," + y + "," + z + "]\n");
		}
		int ix = Nx - 1;
		while (xAxis[ix] > x) {
			ix--;
		}
		int iy = Ny - 1;
		while (yAxis[iy] > y) {
			iy--;
		}
		int iz = Nz - 1;
		while (zAxis[iz] > z) {
			iz--;
		}
		double xd = (x - xAxis[ix]) / (xAxis[ix + 1] - xAxis[ix]);
		double yd = (y - yAxis[iy]) / (yAxis[iy + 1] - yAxis[iy]);
		double zd = (z - zAxis[iz]) / (zAxis[iz + 1] - zAxis[iz]);
		double c00 = grid[ix][iy][iz] * (1 - xd) + grid[ix + 1][iy][iz] * xd;
		double c10 = grid[ix][iy + 1][iz] * (1 - xd) + grid[ix + 1][iy + 1][iz] * xd;
		double c01 = grid[ix][iy][iz + 1] * (1 - xd) + grid[ix + 1][iy][iz + 1] * xd;
		double c11 = grid[ix][iy + 1][iz + 1] * (1 - xd) + grid[ix + 1][iy + 1][iz + 1] * xd;
		double c0 = c00 * (1 - yd) + c10 * yd;
		double c1 = c01 * (1 - yd) + c11 * yd;
		double c = c0 * (1 - zd) + c1 * zd;
		return c;
	}


	private void readMRC(String MRC_filename) {
		try {
			DataInputStream inputStream = new DataInputStream(new FileInputStream(new File(MRC_filename)));
			Nx = myReadInt(inputStream);
			Ny = myReadInt(inputStream);
			Nz = myReadInt(inputStream);
			dataType = myReadInt(inputStream);
			startX = myReadInt(inputStream);
			startY = myReadInt(inputStream);
			startZ = myReadInt(inputStream);
			intervalX = myReadInt(inputStream);
			intervalY = myReadInt(inputStream);
			intervalZ = myReadInt(inputStream);
			dimX = myReadFloat(inputStream);
			dimY = myReadFloat(inputStream);
			dimZ = myReadFloat(inputStream);
			angX = myReadFloat(inputStream);
			angY = myReadFloat(inputStream);
			angZ = myReadFloat(inputStream);
			if (myReadInt(inputStream) != 1) {
				throw new RuntimeException("Should be 1");
			}
			if (myReadInt(inputStream) != 2) {
				throw new RuntimeException("Should be 2");
			}
			if (myReadInt(inputStream) != 3) {
				throw new RuntimeException("Should be 3");
			}
			Dmin = myReadFloat(inputStream);
			Dmax = myReadFloat(inputStream);
			Dmean = myReadFloat(inputStream);
			mySkip(inputStream, 108);
			Ox = myReadFloat(inputStream);
			Oy = myReadFloat(inputStream);
			Oz = myReadFloat(inputStream);
			mySkip(inputStream, 816);
			grid = new float[Nx][Ny][Nz];
			for (int iz = 0; iz < Nz; iz++) {
				for (int iy = 0; iy < Ny; iy++) {
					for (int ix = 0; ix < Nx; ix++) {
						grid[ix][iy][iz] = myReadFloat(inputStream);
					}
				}
			}
			inputStream.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	private void mySkip(DataInputStream inputStream, int NtoSkip) throws IOException {
		for (int c = 0; c < NtoSkip; c++) {
			inputStream.readByte();
		}
	}

	private float myReadFloat(DataInputStream inputStream) throws IOException {
		return Float.intBitsToFloat(myReadInt(inputStream));
	}

	private int myReadInt(DataInputStream inputStream) throws IOException {
		byte a = inputStream.readByte();
		byte b = inputStream.readByte();
		byte c = inputStream.readByte();
		byte d = inputStream.readByte();
		return (((d & 0xff) << 24) | ((c & 0xff) << 16) | ((b & 0xff) << 8) | (a & 0xff)); // redundant var 'v'.
	}

	/**
	 * return actual MRC coordinates from normalized x,y,z coords.
	 * accept int coords normalized for 'grid'
	 * @return array of floats.
	 */
	public float[] getActualCoords(int x, int y, int z) {
		return new float[] {xAxis[x],yAxis[y],zAxis[z]};
	}


}
