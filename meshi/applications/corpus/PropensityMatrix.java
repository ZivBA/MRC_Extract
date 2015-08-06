package meshi.applications.corpus;

import java.util.StringTokenizer;

import meshi.util.Command;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.file.File2StringArray;

public class PropensityMatrix implements KeyWords{
	
	private final double  TORSION_TICKS = (5.0 * (Math.PI/180.0));
	private final double  BRING_ARRAY_TO_ZERO = (Math.PI/TORSION_TICKS);
	private double[][][] prop;
	private double[][] prePro;
	
	public PropensityMatrix(CommandList commands, String fileName) {
	    Command command = commands.firstWord(PARAMETERS_DIRECTORY);
	    String parametersDirectory = command.secondWord(); 
		ReadPropensityFile(parametersDirectory+"/"+fileName);
	}
	
	
	public double propVal(int res, double phi, double psi) {		
		return prop[res][(int) Math.round(phi/TORSION_TICKS+BRING_ARRAY_TO_ZERO)][(int) Math.round(psi/TORSION_TICKS+BRING_ARRAY_TO_ZERO)];
	}
	
	public double preProVal(double phi, double psi) {		
		return prePro[(int) Math.round(phi/TORSION_TICKS+BRING_ARRAY_TO_ZERO)][(int) Math.round(psi/TORSION_TICKS+BRING_ARRAY_TO_ZERO)];
	}
	
	private void ReadPropensityFile(String fileName) {
		int numberOfTicks = (int) Math.round(2*Math.PI/TORSION_TICKS)+1;
		prop = new double[20][numberOfTicks][numberOfTicks];
		prePro = new double[numberOfTicks][numberOfTicks];
		String[] strings = File2StringArray.f2a(fileName);
		int lineCounter = 0;
		int phiC = 0;
		for (int res = 0; res<20 ; lineCounter++) {
			if (!strings[lineCounter].startsWith("#")) {
				int psiC = 0;
				StringTokenizer st = new StringTokenizer(strings[lineCounter]);
				while (st.hasMoreTokens()) {
					prop[res][phiC][psiC] = Double.valueOf(st.nextToken());
					psiC++;
				}
				phiC++;
				if (phiC==numberOfTicks) {
					phiC = 0;
					res++;
				}
			}
		}
		for ( ; lineCounter<strings.length ; lineCounter++) {
			if (!strings[lineCounter].startsWith("#")) {
				int psiC = 0;
				StringTokenizer st = new StringTokenizer(strings[lineCounter]);
				while (st.hasMoreTokens()) {
					prePro[phiC][psiC] = Double.valueOf(st.nextToken());
					psiC++;
				}
				phiC++;
				if (phiC==numberOfTicks) {
					break;
				}
			}
		}
	}

}
