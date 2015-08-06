package meshi.util.dssp;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Vector;

/**
 * Working with the file format that describe the crystal contacts from the CryCo site:
 * http://ligin.weizmann.ac.il/~lpgerzon/cryco5.0/cryco
 * 
 * @author Nir Kalisman
 */

public class CrystalContacts {
	
	private Vector<Integer> contacts = new Vector<Integer>(); 
	
	/**
	 * Constructing the data from the file "fileName" which must be a CryCo file.
	 * All non-solvent contacts under "contactMaxDis" are considered. 
	 *
	 * @param fileName
	 * @param contactMaxDis
	 */
	
	public CrystalContacts(String fileName, double contactMaxDis) {
		String line = "";
		try {
		BufferedReader br = new BufferedReader(new FileReader(fileName));
	    line = br.readLine(); // Skipping the first header line
	    line = br.readLine();
	    while (line != null) {
	    	if (!line.substring(1,4).trim().equalsIgnoreCase("HOH")) { // We don't care about water
	    		Integer resNum = new Integer(line.substring(4,9).trim());
		    	if (!line.substring(30,33).trim().equalsIgnoreCase("HOH")) { // We don't care about water
		    		double dis = (new Double(line.substring(47,51).trim())).doubleValue();
		    		if (dis<contactMaxDis)
		    			if (!isContact(resNum))
		    				contacts.add(resNum);		    		
		    	}
	    	}
	    	line = br.readLine();
	    }
	    br.close();
		}
		catch (Exception e) {
			System.out.println("an error reading line:\n"+line);
			throw new RuntimeException(e);
		}
	}
	
	/**
	 * Return the residue numbers of contacting residues.
	 */
	public int[] getContacts() {
		int[] results = new int[contacts.size()];
		for (int c=0; c<results.length ; c++)
			results[c] = contacts.get(c).intValue();
		return results;
	}
	
	
	/**
	 * If residue 'resNum' is in contact, returns true. 
	 */
	public boolean isContact(int resNum) {
		for (int c=0; c<contacts.size() ; c++)
			if (contacts.get(c).intValue()==resNum)
				return true;
		return false;
	}
	
}
