package meshi.util.crossLinking;

import meshi.util.file.File2StringArray;

/**
 * Simple reading of a sequence from a FASTA file. 
 * @author Nir 
 **/
public class FastaSeq {

    private String fileName = null;
    private String title = null;
	private String seq = null;
	
	public FastaSeq(String fileName) {
    	this.fileName = fileName;
        String[] lines = File2StringArray.f2a(fileName);
        String title = null;
        String seq = "";
        for (int c=0 ; c<lines.length ; c++) {
        	if (!lines[c].startsWith("#")) {  // Skipping comments
        		if (lines[c].startsWith(">")) { // found a title
        			title = lines[c].substring(1).trim();
        		} 
        		else {
        			seq += lines[c].trim();
        		}
        	}
        }
		this.title = title;
		this.seq = seq;
	}
	
	public String title() {
		return title;
	}

	public String seq() {
		return seq;
	}

	public String fileName() {
		return fileName;
	}	
}
