package meshi.util.Dali;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.Vector;

import meshi.applications.loopBuilding.applications.CompleteMissingResidues;
import meshi.applications.loopBuilding.applications.SimpleHomologyModeling;
import meshi.applications.minimizeProtein.ExtendedAtomsProtein;
import meshi.molecularElements.Atom;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.MinimizerException;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.file.File2StringArray;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;

public class PairwiseLiteDali implements Residues, Filter {
	
	/**
	 * This class reads and parses the structural alignment data that 'Pairwise' DaliLight V3.1
	 * produces for the first hit.  
	 * 
	 * @param filename
	 */
	
	protected String rangeFileName = null;
	protected double Zscore = -999;	
	protected int alignmentLength = -999;
	protected double seqIdentity = -999;
	protected int[][] queryRanges = null;
	protected int[][] templateRanges = null;
	protected String queryString = null;
	protected String templateString = null;
	private int[] queryMatch =null;
	private boolean[] okResidues = null;
	
	public PairwiseLiteDali(String filename) {
		rangeFileName = filename;
		System.out.println("Reading the top hit in LiteDali file: " + rangeFileName);
		Vector<String> ranges = getTextFromDaliRangeFile(filename , 1);
		parseHeaderLine(ranges.firstElement());
		queryRanges = new int[ranges.size()-1][2];
		templateRanges = new int[ranges.size()-1][2];
		for (int c=0; c<ranges.size()-1 ; c++) {
			StringTokenizer st = new StringTokenizer(ranges.get(c+1));
			st.nextToken();st.nextToken();st.nextToken();
			st.nextToken();st.nextToken();st.nextToken();
			st.nextToken();st.nextToken();st.nextToken();
			st.nextToken();st.nextToken();
			queryRanges[c][0] = Integer.parseInt(st.nextToken());
			st.nextToken();st.nextToken();
			queryRanges[c][1] = Integer.parseInt(st.nextToken());
			st.nextToken();st.nextToken();
			templateRanges[c][0] = Integer.parseInt(st.nextToken());
			st.nextToken();st.nextToken();
			templateRanges[c][1] = Integer.parseInt(st.nextToken());			
			System.out.println("["+queryRanges[c][0]+","+queryRanges[c][1]+"] - " + 
					"["+templateRanges[c][0]+","+templateRanges[c][1]+"]");
		}
	}
	
	public String makeAlignmentString(String fullQuerySeq, Protein nat, Protein template, int queryStartRes) {
		queryMatch = new int[queryStartRes+fullQuerySeq.length()+9999/*Just to pad it in the end*/];
		for (int c=0 ; c<queryMatch.length ; c++) {
			queryMatch[c] = -999;
		}
		
		// Matching the ranges
		for (int c=0 ; c<queryRanges.length ; c++) {
			int resTemplate = templateRanges[c][0];
			for (int res=queryRanges[c][0] ; res<=queryRanges[c][1] ; res++) {
				if ((nat.residue(res).ca()!=null) && (template.residue(resTemplate).ca()!=null)) {
					queryMatch[res] = resTemplate;
					resTemplate++;
				}
				else if ((nat.residue(res).ca()==null) && (template.residue(resTemplate).ca()==null)) {
					resTemplate++;
				}
				else if ((nat.residue(res).ca()==null) && (template.residue(resTemplate).ca()!=null)) {
				}
				else if ((nat.residue(res).ca()!=null) && (template.residue(resTemplate).ca()==null)) {
					resTemplate++;
					res--;
				}
			}
		}
		
		// Matching between the ranges
		for (int c=0 ; c<(queryRanges.length-1) ; c++) {
			int gapLength = Math.min(queryRanges[c+1][0]-queryRanges[c][1]-1 , 
					templateRanges[c+1][0]-templateRanges[c][1]-1);
			// first from the Nterm
			for (int res=1 ; res<=(gapLength/2) ; res++) {
				queryMatch[queryRanges[c][1]+res] = templateRanges[c][1] + res;				
			}
			// then from the Cterm
			for (int res=gapLength ; res>(gapLength/2) ; res--) {
				queryMatch[queryRanges[c+1][0]-(gapLength-res+1)] = templateRanges[c+1][0]-(gapLength-res+1);				
			}
		}
		
		// Extending from the Nterm
		for (int res=queryRanges[0][0]-1 ; res>=0 ; res--) {
			queryMatch[res] = templateRanges[0][0] - (queryRanges[0][0]-res);
		}
		
		// Extending from the Cterm
		for (int res=queryRanges[queryRanges.length-1][1]+1 ; res<queryMatch.length ; res++) {
			queryMatch[res] = templateRanges[queryRanges.length-1][1] + (res-queryRanges[queryRanges.length-1][1]);
		}
		
		// Making the strings
		String query = "";
		String templ = "";
		String match = "";
		for (int c=0 ; c<fullQuerySeq.length() ; c++) {
			int res = queryStartRes+c;
			// A gap in template
			if (queryMatch[res]==-999) {
				query += fullQuerySeq.charAt(c);
				templ += "-";
			}
			else { 
				query += fullQuerySeq.charAt(c);
				templ += getResTypeFromTemplate(template, queryMatch[res]);
				if (queryMatch[res+1]-queryMatch[res]!=1) { // Put a gap in the query
					for (int cc=1 ; cc<queryMatch[res+1]-queryMatch[res] ; cc++) {
						query += "-";
						templ += getResTypeFromTemplate(template, queryMatch[res]+cc);
					}
				}
			}
		}
		// Makeing the match string
		for (int c=0; c<query.length() ; c++) {
			if (query.charAt(c)==templ.charAt(c)) {
				match += query.charAt(c);
			}
			else {
				match += " ";
			}
		}
		
		String output = queryStartRes + "\n" +
		Math.max(template.firstResidue() , queryMatch[queryStartRes]) + "\n" +
		query + "\n" + 
		match + "\n" + 
		templ + "\n"; 
		
		queryString = query;
		templateString = templ;
		
		// updating 'queryMatch' and 'queryMatch' for further use
		for (int c=0; c<queryMatch.length ; c++) {
			if (queryMatch[c]!=-999) {
				if (getResTypeFromTemplate(template, queryMatch[c]).equals("X")) {
					queryMatch[c]=-999;
				}
			}
		}
		okResidues = new boolean[queryMatch.length];
		
		return output;
	}
	
	
	private Vector<String> getTextFromDaliRangeFile(String filename ,int  daliHit) {
		String line;
		Vector<String> returnStrings = new Vector<String>();
		String lookForString = "Alignment number = " + daliHit;
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			line = br.readLine(); 
			while (line.indexOf(lookForString)==-1) {
				line = br.readLine();
			}
			returnStrings.add(line);
			line = br.readLine();
			while ((line!=null) && (line.indexOf("Alignment number =")==-1) && (line.indexOf("<=>")>-1)) {
				returnStrings.add(line);
				line = br.readLine();
			}
			br.close();
		}
		catch (Exception e) {
			System.out.println("\n\nAn error reading file: "+filename + 
					"\nThis is not a \'range.txt\' file produced by DaliLite, or it does not contain a line with \'alignment number\'\n\n");
			throw new RuntimeException(e);
		}
		return returnStrings;		
	}

	private void parseHeaderLine(String header) {
		StringTokenizer st = new StringTokenizer(header);
		st.nextToken();st.nextToken();st.nextToken();st.nextToken();
		st.nextToken();st.nextToken();st.nextToken();
		String tmpString = st.nextToken();
		String ZscoreString = tmpString.substring(0, tmpString.length()-1); // removing the end comma
		Zscore = Double.parseDouble(ZscoreString);
		System.out.println("Zscore: " + Zscore);
		st.nextToken();st.nextToken();st.nextToken();st.nextToken();
		tmpString = st.nextToken();
		String alignmentLengthString = tmpString.substring(0, tmpString.length()-1); // removing the end comma
		alignmentLength = Integer.parseInt(alignmentLengthString);
		System.out.println("Alignment Length: " + alignmentLength);
		st.nextToken();st.nextToken();st.nextToken();
		seqIdentity = Double.parseDouble(st.nextToken());
		System.out.println("Sequence identity: " + seqIdentity + "%");
	}
    
    
    /**
     * This method return the residue type (in one letter) of 'resNum' in 'prot'. If This residue is not in 
     * prot it returns "X" if it is in the middle of the protein or "-" if it is past its termini.
     */
    public static String getResTypeFromTemplate(Protein prot, int resNum) {
    	if ((prot.residue(resNum)!=null) && (prot.residue(resNum).ca()!=null)) {
    		return Residue.nameOneLetter(prot.residue(resNum).type);
    	}
    	int firstResNum = prot.firstResidue();
    	int lastResNum = prot.lastResidue();
    	if ((resNum<=lastResNum) && (resNum>=firstResNum)) {
    		return "X";
    	}
    	else {
    		return "-";
    	}    	
    }

    // Run like this: java PairwiseLiteDali ProteinName SEQ_file Native_file Template_file Dali_file 
    public static void main(String[] args) throws MinimizerException, LineSearchException, IOException {
    	String path = "/home/nirka/projects/directRef/"; // if not empty string must ends with '/'
    	String commandFileName = "commands";

    	String protID = args[0].trim();
    	MeshiProgram.initRandom(777);
    	Protein nat = new Protein(path + "Natives/"+protID+".pdb" , new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
    	Protein templ = new Protein(path + "Templates/Template_"+protID+".pdb"  , new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
    	PairwiseLiteDali align = new PairwiseLiteDali(path + "DALI/"+protID+"/ranges.txt");

    	// Getting the sequence:
    	String[] fastaSeq = File2StringArray.f2a(path + "Sequences/"+protID.substring(0,5)+".seq" );
    	String fullSeq = "";
    	for (int c=1 ; c<fastaSeq.length ; c++) 
    		fullSeq += fastaSeq[c];
    	int firstResNum = Math.max(1, nat.firstResidue()-3);
    	int lastResNum = Math.min(fullSeq.length(), nat.lastResidue()+3);
    	String seq = fullSeq.substring(firstResNum-1, lastResNum);
    	
    	
    	// making the alignment string
    	String alignmentString = path + "Templates/Template_"+protID+".pdb" + "\n" +
    	align.makeAlignmentString(seq,nat,templ,firstResNum) + 
    	path + "PreModels/"+protID+".pdb" + "\n\n";
    	System.out.println("\n\nAlignment:\n----------\n"+alignmentString);
    	
    	// Writing it to disk
    	try{
    		BufferedWriter bw = new BufferedWriter(new FileWriter( path + "Alignments/"+protID+".align" ));
    		bw.write(alignmentString);
    	    bw.close();
    	}
    	catch(Exception e) {
    	    throw new RuntimeException(e.getMessage());
    	} 
    	
    	// Making the segmod file
    	String segmodString = align.getQueryString().length() + " " + "Template_" + protID + " XXX\n";
    	for (int qCounter=0; qCounter<align.getQueryString().length() ; qCounter++) {
    		segmodString += " " + align.getTemplateString().charAt(qCounter);
    		if ((qCounter % 30) == 29) {
    			segmodString += "\n";
    		}
    	}
    	if ((align.getQueryString().length() % 30) != 0) {
    		segmodString += "\n";
    	}
    	segmodString += align.getQueryString().length() + " " + "SEGMOD_" + protID + " YYY\n";
    	for (int qCounter=0; qCounter<align.getQueryString().length() ; qCounter++) {
    		segmodString += " " + align.getQueryString().charAt(qCounter);
    		if ((qCounter % 30) == 29) {
    			segmodString += "\n";
    		}
    	}
    	if ((align.getQueryString().length() % 30) != 0) {
    		segmodString += "\n";
    	}
    	try{
    		BufferedWriter bw = new BufferedWriter(new FileWriter( path + "PreModels/"+protID+".ali" ));
    		bw.write(segmodString);
    	    bw.close();
    	}
    	catch(Exception e) {
    	    throw new RuntimeException(e.getMessage());
    	} 
    	
    	
    	// Homology modeling
    	String[] argsString_homology = {commandFileName , path + "Alignments/"+protID+".align" };
    	SimpleHomologyModeling.main(argsString_homology);
    	
    	// Adding missing parts, but not if the gap is more than 10 residues
    	Protein query = CompleteMissingResidues.addResiduesToProtein(new ExtendedAtomsProtein(path + "PreModels/"+protID+".pdb",DO_NOT_ADD_ATOMS) , 
				firstResNum, seq, firstResNum, firstResNum+seq.length()-1);
    	int maxGap = 10;
    	int includeTerm = 2;
    	for(int c=0; c<align.getOkResidues().length ; c++) {
    		align.getOkResidues()[c] = false;
    		if (align.getQueryMatch()[c]!=-999) {
    			align.getOkResidues()[c]=true;
    		}
    		else {
    			// search left
    			int gapLeftLength = 0;
    			while ((c-gapLeftLength>-1) && (align.getQueryMatch()[c-gapLeftLength-1]==-999)) {
    				gapLeftLength++;
    			}
    			// search right
    			int gapRightLength = 0;
    			while ((c+gapRightLength<align.getOkResidues().length) && (align.getQueryMatch()[c+gapRightLength]==-999)) {
    				gapRightLength++;
    			}
    			if ((gapRightLength+gapLeftLength<maxGap) || (gapLeftLength<includeTerm) ||
    					(gapRightLength<=includeTerm)) {
    				align.getOkResidues()[c]=true;
    			}
    		}
    	}
    	Protein filteredQuery = new Protein(query.atoms().filter(align) ,
    			new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
    	filteredQuery.atoms().print(new MeshiWriter(path + "PreModels/"+protID+"_full.pdb"));
    	
//    	// Beautifying before loop run:
//    	RefinePreModel.refine(new CommandList(commandFileName), path + "PreModels/"+protID+"_full.pdb" ,
//    			path + "PreModels/"+protID+".pdb" ,
//    			path + "Models/"+protID+".pdb");
    	
    }

    
    protected int[] getQueryMatch() {return queryMatch;}
    protected boolean[] getOkResidues() {return okResidues;}
    protected String getQueryString() {return queryString;}
    protected String getTemplateString() {return templateString;}
    
    public boolean accept(Object obj) {
    	Atom atom = (Atom) obj;
    	return okResidues[atom.residueNumber()];
    }
    
}
