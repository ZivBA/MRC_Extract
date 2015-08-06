package meshi.energy.ROT1solvation.parameters;

import java.io.File;
import java.util.StringTokenizer;

import meshi.geometry.Distance;
import meshi.util.file.File2StringArray;

/**
 * 
 * The purpose of implementations for this class are to give the solvation value array through the
 * 'getROT1SolvValueArray' method. Each row in the array is for a certain atomType. 
 * For each 'atomType' the row is an array of the solvation values for cnc's of: 
 * {0,1,2,3,.....}, i.e. equal spacing of 1 for the cnc values.
 * Solvation for (cnc<0) or (cnc>=cncMax) is 0.0
 * 
 * Implementaions should also say (through 'filterDisForRelevance') if a certain distance is relevant 
 * for the CNC calculation. 
 * 
 * Certain sigmoid values should also be given through the abstract methods.
 * 
 * @author Nir
 *
 */

public abstract class AbstractROT1Parameters {

	protected double[][] solvEne = null; 
	
	public AbstractROT1Parameters(String parameterFileName) {
		readSolvEne(parameterFileName);
	}
		
	public double[][] getROT1SolvValueArray() {
		return solvEne;
	}

	public double end() {
		return cutoff()+2*sigmoidRange();
	}

	public double p1() {
		return cutoff()-sigmoidRange();
	}

	public double p2() {
		return cutoff()+sigmoidRange();
	}

	public double valAtp1() {
		return 0.99;
	}

	public double valAtp2() {
		return 0.01;
	}

	public boolean filterDisForRelevance(Distance dis) {
		if (dis.atom1().isHydrogen || dis.atom2().isHydrogen)
			return false;
		if (dis.atom1().isBackbone&&dis.atom2().isBackbone)
			return false;
		if ((dis.atom1().residueNumber()>(dis.atom2().residueNumber()-2)) && 
				(dis.atom1().residueNumber()<(dis.atom2().residueNumber()+2)))
			return false;
		return true;
	}
	
	protected void readSolvEne(String filename) {
		if ((new File(filename)).exists()) 
			System.out.println("Reading Parameter file: " + filename);
		else
			throw new RuntimeException("\n\nParameter file not found: " + filename + "\n\n");
		String[] strArray = File2StringArray.f2a(filename);
		solvEne = new double[strArray.length][];
		for (int c=0 ; c<strArray.length ; c++) {
			StringTokenizer st = new StringTokenizer(strArray[c]);
			solvEne[c] = new double[st.countTokens()];
			for (int d=0; st.hasMoreTokens() ; d++)
				solvEne[c][d] = (new Double(st.nextToken())).doubleValue();
		}
	}

	protected abstract double cutoff();
	
	protected abstract double sigmoidRange();

}
