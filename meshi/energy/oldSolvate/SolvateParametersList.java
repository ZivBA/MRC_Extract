package meshi.energy.oldSolvate;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.energy.ParametersList;
import meshi.molecularElements.Atom;
import meshi.util.mathTools.Spline1D;

/**
 * The parameter class required for the SolvateEnergy class. This class reads the parameters from 12 different files.
 * The 12 file names are given as a String array to the constructor. 
 * 
 * The 10 last files, are files that read parameters for the different sigmoid functions used in the calculations 
 * of the carbon and HB indices in SolvateEnergy. See the documentation of the Sigma class to see what are the 
 * parameters required to define a sigmoid. The formats of all these files are similar. Each file contains a
 * 14x14 matrix of doubles (14 is the number of types in Tsai 99'). Value i,j in the matrix, is the relevent sigmoid 
 * property of atom type j on atom type i. For example, the sixth value in the second row, corresponds to a sigmoid 
 * of hydroxyl oxygen (type 6) on a backbone nitrogen (type 2).     
 *
 * The 12 file types:
 * ------------------
 * 1) SolvateMESHI2Tsai.dat - In MESHI there are 190 defined atom types. In the solvate energy, we sometime use 
 * a reduce representation of 14 atom types as defined in Tsai et. al. 99'. The format of this file is:
 * {MESHI atom type (string)} {Tsai type number (int)}
 * {MESHI atom type (string)} {Tsai type number (int)}  
 * {MESHI atom type (string)} {Tsai type number (int)}
 *  ...
 * {MESHI atom type (string)} {Tsai type number (int)}  
 *
 * 2) SolvateExtResParameters.dat - 1D spline parameters for the conversion of the environment indices into 
 * energies. A spline is created for each of the of the 190 atom types. The format of this file is:
 * {MESHI atom type (string)} {1D spline initialization string (see Spline1D.java)}
 * {MESHI atom type (string)} {1D spline initialization string (see Spline1D.java)}
 * {MESHI atom type (string)} {1D spline initialization string (see Spline1D.java)}
 *  ...
 * {MESHI atom type (string)} {1D spline initialization string (see Spline1D.java)}
 *
 * 3) SolvateExtResCend.dat - A 14x14 matrix of doubles. These values are the 'end' values in the sigmoid that
 * calculates the carbon index. 
 *
 * 4) SolvateExtResCp1.dat - A 14x14 matrix of doubles. These values are the 'p1' values in the sigmoid that
 * calculates the carbon index. 
 *
 * 5) SolvateExtResCp2.dat - A 14x14 matrix of doubles. These values are the 'p2' values in the sigmoid that
 * calculates the carbon index. 
 *
 * 6) SolvateExtResCvalAtp1.dat - A 14x14 matrix of doubles. These values are the 'valAtp1' values in the sigmoid that
 * calculates the carbon index. 
 *
 * 7) SolvateExtResCvalAtp2.dat - A 14x14 matrix of doubles. These values are the 'valAtp2' values in the sigmoid that
 * calculates the carbon index. 
 *
 * 8) SolvateExtResHBend.dat - A 14x14 matrix of doubles. These values are the 'end' values in the sigmoid that
 * calculates the HB index. 
 *
 * 9) SolvateExtResHBp1.dat - A 14x14 matrix of doubles. These values are the 'p1' values in the sigmoid that
 * calculates the HB index. 
 *
 * 10) SolvateExtResHBp2.dat - A 14x14 matrix of doubles. These values are the 'p2' values in the sigmoid that
 * calculates the HB index. 
 *
 * 11) SolvateExtResHBvalAtp1.dat - A 14x14 matrix of doubles. These values are the 'valAtp1' values in the sigmoid that
 * calculates the HB index. 
 *
 * 12) SolvateExtResHBvalAtp2.dat - A 14x14 matrix of doubles. These values are the 'valAtp2' values in the sigmoid that
 * calculates the HB index. 
 *
 **/

public class SolvateParametersList extends ParametersList {
	
	public final int NTsai = 14; // The number of atom types used in Tsai 99'. Any Hydrogen is type 14.
    public final int[] atomicTypeConverter; // From MESHI atom types to Tsai 99'.
    public final double maxEnd; // The maximal distance where any sigmoid is not zero. This value must be less than the Rmax in the distance matrix. 
    public final Spline1D[] atomTypeSplines;
    public final double[][] Cend;     
    public final double[][] Cp1;
    public final double[][] Cp2;
    public final double[][] CvalAtp1;
    public final double[][] CvalAtp2;
    public final double[][] HBend;     
    public final double[][] HBp1;
    public final double[][] HBp2;
    public final double[][] HBvalAtp1;
    public final double[][] HBvalAtp2;

/**
 * The parameter to the constructor is an array of 12 Strings, giving the 12 files required as parameters.
 * See the MeshiPotential class for a detailed list.
 **/    
    public SolvateParametersList(String[] parameterFiles) {
	super();
	int maxAtomType=-1;
	int tmp;
	BufferedReader br;
    StringTokenizer stok;
    String line = "";
	// Converting the 190 MESHI atom types to the 14 mentioned in Tsai 99'
	System.out.println("Reading solvation parameter file: " + parameterFiles[0]);
	try{
		// first pass on the file - to find the maximal atom type
		br = new BufferedReader(new FileReader(parameterFiles[0]));
	    line = br.readLine();
	    while (line != null) {
	    	stok = new StringTokenizer(line);
	    	tmp = Atom.type(stok.nextToken().trim());
	    	if (tmp>maxAtomType)
	    	   maxAtomType = tmp;
	    	line = br.readLine();
	    }
	    br.close();
	    atomicTypeConverter = new int[maxAtomType+1];
	    for(int c=0 ; c<atomicTypeConverter.length; c++)
	    	atomicTypeConverter[c] = -1; 
		// second pass on the file - reading the new types
		br = new BufferedReader(new FileReader(parameterFiles[0]));
	    line = br.readLine();
	    while (line != null) {
	    	stok = new StringTokenizer(line);
	    	tmp = Atom.type(stok.nextToken().trim());
	    	atomicTypeConverter[tmp] = Integer.valueOf(stok.nextToken().trim()).intValue()-1;
	    	line = br.readLine();
	    }
	    br.close();
	}
	catch(Exception e) {
	    throw new RuntimeException(e.getMessage());
	}
	
	// Reading the spline parameters for the 190 atom types.
	atomTypeSplines = new Spline1D[maxAtomType+1];
	System.out.println("Reading solvation parameter file: " + parameterFiles[1]);
	try{
		br = new BufferedReader(new FileReader(parameterFiles[1]));
	    line = br.readLine();
	    while (line != null) {
	    	stok = new StringTokenizer(line);
	    	tmp = Atom.type(stok.nextToken().trim());
	    	atomTypeSplines[tmp] = new Spline1D(stok);
	    	line = br.readLine();
	    }
	    br.close();
	}
	catch(Exception e) {
	    throw new RuntimeException(e.getMessage());
	}
		
	Cend = new double[NTsai][NTsai];
	readSigmoidValueFile(Cend,parameterFiles[2]);

	Cp1 = new double[NTsai][NTsai];
	readSigmoidValueFile(Cp1,parameterFiles[3]);

	Cp2 = new double[NTsai][NTsai];
	readSigmoidValueFile(Cp2,parameterFiles[4]);

	CvalAtp1 = new double[NTsai][NTsai];
	readSigmoidValueFile(CvalAtp1,parameterFiles[5]);

	CvalAtp2 = new double[NTsai][NTsai];
	readSigmoidValueFile(CvalAtp2,parameterFiles[6]);

	HBend = new double[NTsai][NTsai];
	readSigmoidValueFile(HBend,parameterFiles[7]);

	HBp1 = new double[NTsai][NTsai];
	readSigmoidValueFile(HBp1,parameterFiles[8]);

	HBp2 = new double[NTsai][NTsai];
	readSigmoidValueFile(HBp2,parameterFiles[9]);

	HBvalAtp1 = new double[NTsai][NTsai];
	readSigmoidValueFile(HBvalAtp1,parameterFiles[10]);

	HBvalAtp2 = new double[NTsai][NTsai];
	readSigmoidValueFile(HBvalAtp2,parameterFiles[11]);
		
	// Finding the largest ending distance
	double maxEndTmp = -1.0;
	for (int c1=0 ; c1<NTsai ; c1++) 
	for (int c2=0 ; c2<NTsai ; c2++) {
	   if (Cend[c1][c2]>maxEndTmp)
	      maxEndTmp = Cend[c1][c2];
	   if (HBend[c1][c2]>maxEndTmp)
	      maxEndTmp = HBend[c1][c2];
	}
	maxEnd = maxEndTmp;

    }
    
    public Parameters createParameters(String s) {
    	throw new RuntimeException("This method should not be called");
    }
    
    public Parameters parameters(Object obj) {
    	    	throw new RuntimeException("This method should not be called");
    }
    
    
    // Reading a 14x14 file into the relevent parameter.
	private void readSigmoidValueFile(double[][] ar, String filename) {
		BufferedReader br;
		StringTokenizer stok;
		String line = "";
		System.out.println("Reading solvation parameter file: " + filename);
		try{
			br = new BufferedReader(new FileReader(filename));
			for (int ind1=0 ; ind1<NTsai ; ind1++) {
				line = br.readLine();
				if (line==null) 
					throw new RuntimeException("In " + filename + " there should be a " + NTsai + "x" +
					NTsai + " array of doubles");
				stok = new StringTokenizer(line);
				for (int ind2=0 ; ind2<NTsai ; ind2++) 
					ar[ind1][ind2] = Double.valueOf(stok.nextToken().trim()).doubleValue();
			}
			br.close();
		}
		catch(Exception e) {
			throw new RuntimeException(e.getMessage());
		}
	}    
}
