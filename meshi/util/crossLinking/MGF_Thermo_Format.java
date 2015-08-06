package meshi.util.crossLinking;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.StringTokenizer;
import java.util.Vector;

public class MGF_Thermo_Format {

	private IonVector iVec = null;
	 
	public MGF_Thermo_Format(String fileName) {
    	double mass = -1;
    	double intensity = -1;
    	int chargeState = -1;
    	int scan = -1;
    	double retentionTime = -1;
    	String ionFileName = null;
    	boolean readingIonList = false;
    	Vector<double[]> ms2SpectrumCollection = new Vector<double[]>();
    	iVec = new IonVector();
    	String line="";
    	try {
    		BufferedReader fMGF = new BufferedReader(new FileReader(fileName));
    		line = fMGF.readLine();
    		while (line != null) {
    			if (line.startsWith("END IONS")) {
    				// preparing the ms2 spectrum.
    				double[][] ms2Spectrum = new double[ms2SpectrumCollection.size()][2];
    				for (int c=0 ; c<ms2SpectrumCollection.size() ; c++) {
    					ms2Spectrum[c][0] = ms2SpectrumCollection.get(c)[0];
    					ms2Spectrum[c][1] = ms2SpectrumCollection.get(c)[1];
    				}
    				Ion ion = new Ion(mass, chargeState, scan, ionFileName, intensity, retentionTime, ms2Spectrum);
    				iVec.add(ion);
    				readingIonList = false;
    			}
    			if (readingIonList) {
    				StringTokenizer st = new StringTokenizer(line);
    				double mz = Double.valueOf(st.nextToken());
    				double intensityMS2 = Double.valueOf(st.nextToken());
    				double[] tmpDoubleVec = {mz , intensityMS2};
    				ms2SpectrumCollection.add(tmpDoubleVec);
    			}
    			else {
        			if (line.startsWith("PEPMASS=")) {
        				StringTokenizer tmpST = new StringTokenizer(line);
        				mass = Double.valueOf(tmpST.nextToken().substring(8));
        				intensity = Double.valueOf(tmpST.nextToken());
        			}
        			if (line.startsWith("CHARGE=")) {
        				chargeState = Integer.valueOf(line.substring(7));
        				readingIonList = true;
        				ms2SpectrumCollection = new Vector<double[]>();
        			}    				
        			if (line.startsWith("RTINSECONDS=")) {
        				retentionTime = Double.valueOf(line.substring(12));
        			}    				
        			if (line.startsWith("TITLE=")) {
        				StringTokenizer tmpST = new StringTokenizer(line);
        				ionFileName = fileName;
        				tmpST.nextToken();
        				tmpST.nextToken();        				
        				scan = Integer.valueOf(tmpST.nextToken().substring(5));
        			}    				        			
    			}
    			line = fMGF.readLine();
    		} 
    		fMGF.close();
    	}
    	catch(Exception e) {
    		System.out.println("GGG\n\n\n" + line + "\n\n\nGGGGG\n");
    		throw new RuntimeException(e);
    	}
	}

	public String toString() {
		String out = "";
		for (int c=0 ; c<iVec.size() ; c++) {
			out += (iVec.get(c).scan() + " " + iVec.get(c).intensityMS1() + "\n"); 
		}
		return out;
	}
	
}
