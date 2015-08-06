package meshi.util.crossLinking;

import java.util.LinkedList;
import java.util.Vector;

import meshi.util.file.File2StringArray;

public class MySequenceList extends Vector<MySequence> {
	
    private String fileName = null;

    public MySequenceList(String fileName) {
    	this.fileName = fileName;
        String[] lines = File2StringArray.f2a(fileName);
        String title = null;
        String seq = "";
        for (int c=0 ; c<lines.length ; c++) {
        	if (!lines[c].startsWith("#")) {  // Skipping comments
        		if (lines[c].startsWith(">")) { // found a title
        			if (title != null) {
        				add(new MySequence(title, seq));
        				title = null;
        			}
        			title = lines[c].substring(1).trim();
        			seq = "";
        		}
        		else { // lengthening the sequence
        			seq += processLine(lines[c]);
        		}
        	}
        }
        // Adding the last sequence
        if (title!=null) {
        	add(new MySequence(title, seq));
        }      
    }
	
	public String processLine(String rawString) {
		return rawString.trim();		
	}
	
    public String fileName() {
		return fileName;
	}


	/**
     * Returning list of all the instances where 'seq' appears in the sequence list.
     * In the returning list: found[X][0] is the sequence index in the sequence list (zero-based), 
     * and found[X][1] is the begining residue in the sequence (1-based).  
     **/
    public int[][] findInSequences(String querySeq) {
    	LinkedList<int[]> found = new LinkedList<int[]>();
    	for (int c=0 ; c<size() ; c++) {
    		MySequence mySeq = get(c);
    		int fromIndex = -1;
    		while (fromIndex<mySeq.seq().length()) {
    			int foundInd = myIndexOf(mySeq.seq(),querySeq, fromIndex);
    			if (foundInd==-1) {
    				fromIndex = mySeq.seq().length() + 1;
    			} 
    			else {
    				int[] tmp = new int[2];
    				tmp[0] = c;
    				tmp[1] = foundInd + 1;    				
    				found.add(tmp);
    				fromIndex = foundInd + 1;
    			}
    		}
    	}
    	int[][] foundAr = new int[found.size()][];
    	for (int c=0 ; c<foundAr.length ; c++) {
    		foundAr[c] = found.get(c);
    	}
    	return foundAr;
    }
    
    /**
     * I had to write my own 'indexOf' to compansate for I/L switches
     */
    private int myIndexOf(String origStr , String origQuery, int fromIndex) {
    	if (fromIndex>=origStr.length()) {
    		return -1;
    	}
    	if (fromIndex<0) {
    		fromIndex=0;
    	}
    	String str = origStr.replace("I", "!").replace("L", "!").replace("SAT", "@@@").replace("TAS", "@@@"); 
    	String query = origQuery.replace("L", "!").replace("I", "!").replace("SAT", "@@@").replace("TAS", "@@@"); 
    	for (int c=fromIndex ; c<=(str.length()-query.length()) ; c++) {
    		boolean globalMatch = true;
    		for (int d=0; globalMatch && (d<query.length()) ; d++) {
    			if (str.charAt(c+d)!=query.charAt(d)) {
    				globalMatch = false;
    			}
    		}
    		if (globalMatch) {
    			return c;
    		}
    	}
    	str = str.replace("N", "#").replace("GG", "#"); 
    	query = query.replace("N", "#").replace("GG", "#"); 
    	for (int c=fromIndex ; c<=(str.length()-query.length()) ; c++) {
    		boolean globalMatch = true;
    		for (int d=0; globalMatch && (d<query.length()) ; d++) {
    			if (str.charAt(c+d)!=query.charAt(d)) {
    				globalMatch = false;
    			}
    		}
    		if (globalMatch) {
    			return c;
    		}
    	}    	
    	return -1;
    }
	
}
