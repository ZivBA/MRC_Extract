package meshi.util.crossLinking;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.StringTokenizer;
import java.util.Vector;

public class MGF_Mass_Matrix_Format {

	private Vector<double[]> ms2 = null;
	private IonVector iVec = null;
	 
	public MGF_Mass_Matrix_Format(String fileName) {
    	ms2 = new Vector<double[]>();
    	double mass = -1;
    	int chargeState = -1;
    	int scan = -1;
    	double retentionTime = -1;
    	String ionFileName = null;
    	boolean readingIonList = false;
    	double sumIons = 0.0;
    	double numIons = 0.0;
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
    				Ion ion = new Ion(mass, chargeState, scan, ionFileName, sumIons, retentionTime, ms2Spectrum);
    				iVec.add(ion);

    				double[] tmp = {mass , chargeState, sumIons , scan};
					ms2.add(tmp);
    				readingIonList = false;
    				sumIons = 0.0;
    				numIons = 0.0;
    			}
    			if (readingIonList) {
    				StringTokenizer st = new StringTokenizer(line);
    				double mz = Double.valueOf(st.nextToken());
    				double intensity = Double.valueOf(st.nextToken());
    				double[] tmpDoubleVec = {mz , intensity};
    				ms2SpectrumCollection.add(tmpDoubleVec);
    				numIons++;
    				sumIons += intensity;
    			}
    			else {
        			if (line.startsWith("PEPMASS=")) {
        				mass = Double.valueOf(line.substring(8));
        				readingIonList = true;
        				ms2SpectrumCollection = new Vector<double[]>();
        			}
        			if (line.startsWith("CHARGE=")) {
        				chargeState = Integer.valueOf(line.substring(7,8));
        			}    				
        			if (line.startsWith("RTINSECONDS=")) {
        				retentionTime = Double.valueOf(line.substring(12));
        			}    				
        			if (line.startsWith("SCANS=")) {
        				scan = Integer.valueOf(line.substring(6));
        			}    	
        			if (line.startsWith("TITLE=File:")) {
        				StringTokenizer tmpST = new StringTokenizer(line);
        				ionFileName = tmpST.nextToken().substring(11);
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
		for (int c=0 ; c<ms2.size() ; c++) {
//			out += (ms2.get(c)[0] + " " + ms2.get(c)[1] + " " + ms2.get(c)[2] + " " + ms2.get(c)[3] + "\n"); 
			out += (ms2.get(c)[3] + " " + ms2.get(c)[2] + "\n"); 
		}
		return out;
	}
	
	public static void main(String[] args) {
//		MGF_Mass_Matrix_Format mgf = new MGF_Mass_Matrix_Format("C:\\Inetpub\\wwwroot\\ISB\\data\\101003_NKalisman_Frydman_NK_20.mgf");
		MGF_Thermo_Format mgf = new MGF_Thermo_Format("C:\\Inetpub\\wwwroot\\ISB\\data\\101003_NKalisman_Frydman_NK_20_thermo.mgf");
		System.out.print(mgf);
	}
	
}
