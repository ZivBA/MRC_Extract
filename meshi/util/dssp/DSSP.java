package meshi.util.dssp;
import java.io.BufferedReader;
import java.io.FileReader;

import meshi.molecularElements.Residue;
import meshi.parameters.Residues;
import meshi.sequences.AccesibilitySequence;
import meshi.sequences.ResidueSequence;
import meshi.sequences.SecondaryStructureSequence;
import meshi.sequences.Sequence;
import meshi.sequences.SequenceList;

public class DSSP implements Residues {
    protected String fileName;
    private char[] aa;
    private int[] resNum;
    private char[] ss;
    private char[] chainID;
    private double[] solvACC;
    private double[] fullSol = {115,135,150,190,210,
				75 ,195,175,200,170,
				185,160,145,180,225,
				115,140,155,255,230};
    private double hbPrecentage;
    private double parallelPrecentage;
    private double antyparallelPrecentage;

    public void printPresentage(){
	System.out.println("Total hb: "+hbPrecentage);
	System.out .println("parallel: "+parallelPrecentage );
	System.out .println("AntiParallel: "+antyparallelPrecentage );

	System.out.println("Number Of Betta Residues: ");
    }

    public DSSP(String dsspFileName)
    {
    	int zvl = ALA; // force the reading of "meshi.parameters.Residues"
    	this.fileName = dsspFileName;
    	readDSSP();
    }
	
    /*
     *This function returns true if the DSSP assignment of residue number 'res' 
     *appears in the given char list. Else false is returned.
     */
    public boolean inList(int res, char chain, char[] SSlist) {
    	char ch;
    	int i;

    	ch = SSofRes(res,chain);
    	for (i=0 ; i<SSlist.length ; i++) {
    		if (ch == SSlist[i])
    			return true;
    	}
    	return false;
    }


    /*
     *This function returns true if the residue given is not in the edge of a structure. 
     *i.e. the residues +/-1 has the same SS as it has.
     */
    public boolean notInEdge(int res, char chain) {
    	char chp1='-',chm1='-',ch='-';
    	int i;
    	for (i=0; (i<resNum.length) ; i++) {
    		if ((resNum[i] == (res-1)) && ((chainID[i] == chain) || (chain == ' ')))
    			chm1 = ss[i];
    		if ((resNum[i] == (res)) && ((chainID[i] == chain) || (chain == ' ')))
    			ch = ss[i];
    		if ((resNum[i] == (res+1)) && ((chainID[i] == chain) || (chain == ' '))) 
    			chp1 = ss[i];
    	}
    	if (chm1=='-')
    		throw new RuntimeException("\n\nCould not find the chain: " + chain + "    or residue: " + res + "\n\n");
    	return (chp1 == ch) && (ch == chm1);
    }	 	

    /*
     *Returns the SS structure of residue number 'res'.
     *If this residue number doesn't exist - 'X' is returned.
     */
    public char SSofRes(int res, char chain) {
    	int i;
    	for (i=0; i<resNum.length ; i++) {
    		if ((resNum[i] == res) && ((chainID[i] == chain) || (chain == ' ')))
    			return ss[i];
    	}
    	return 'X';
    } 


    /*
     *Returns the AA of residue number 'res'.
     *If this residue number doesn't exist - 'X' is returned.
     */
    public char AAofRes(int res, char chain) {
    	int i;
    	for (i=0; i<resNum.length ; i++) {
    		if ((resNum[i] == res) && ((chainID[i] == chain) || (chain == ' ')))
    			return aa[i];
    	}
    	return 'X';
    } 
	
    /*
     *Returns the solvate accessibilty of residue number 'res'.
     *If this residue number doesn't exist - (-1) is returned.
     */
    public double ACCofRes(int res, char chain) {
    	int i;
    	for (i=0; i<resNum.length ; i++) {
    		if ((resNum[i] == res) && ((chainID[i] == chain) || (chain == ' ')))
    			return solvACC[i];
    	}
    	return -1;
    } 


    /*
     *Returns the relative solvate accessibilty of residue number 'res'.
     *If this residue number doesn't exist - (-1) is returned.
     */
    public double relACCofRes(int res, char chain) {
    	int i;
    	for (i=0; i<resNum.length ; i++) {
    		if ((resNum[i] == res) && ((chainID[i] == chain) || (chain == ' ')))
    			return solvACC[i]/fullSol[Residue.type("" + aa[i])];
    	}
    	return -1;
    }

    public String getACC(){
    	String ans = "";
    	for (int i=0; i<resNum.length ; i++)
    		if (resNum[i] > -999)
    			ans = ans+  solvACC[i]+";";
    	return ans;
    }

    public String getRelativeACC(){
    	String ans = "";
    	for (int i=0; i<resNum.length ; i++)
    		if (resNum[i] > -999){
    			if ((int)((solvACC[i]/fullSol[Residue.type("" + aa[i])])+0.5) >= 1)
    				ans = ans+ "A" ;
    			else
    				ans = ans+"B";
    		}
    	return ans;
    }

    public String getAA(){
    	String ans = "";
    	for (int i=0; i<resNum.length ; i++){
    		if (resNum[i] > -999)
    			ans = ans+  aa[i];
    		else
    			throw new RuntimeException("res "+i+" is -999"+"\n"+
    					ans);
    	}
    	return ans;
    }

    public void printSS() {
    	for (int i=0; i<resNum.length ; i++) {
    		if (resNum[i] > -999)
    			System.out.println(ss[i]);
    		//System.out.println(resNum[i]+ "   " + aa[i] + "   " + ss[i] );
    	}
    }

    public String getSSInOneLine() {
    	String ans = "";
    	for (int i=0; i<resNum.length ; i++) {
    		if (resNum[i] > -999)
    			ans = ans+  ss[i];
    		/*if(ss[i] == 'B' | ss[i] =='E')
	      ans = ans+'E';
	      else if(ss[i] == 'H')
	      ans = ans+'H';
	      else if(ss[i] == 'S' | ss[i] == 'G' | ss[i] == 'I' | ss[i] == 'T')
	      ans = ans+'A';
	      else
	      ans = ans+'C';*/

    		//System.out.println(resNum[i]+ "   " + aa[i] + "   " + ss[i] );
    	}
    	return ans;
    }

    public String getSSInOneLineConventional() {
    	String ans = "";
    	for (int i=0; i<resNum.length ; i++) {
    		if (resNum[i] > -999)

    			if(ss[i] == 'B' | ss[i] =='E')
    				ans = ans+'E';
    			else if(ss[i] == 'H')
    				ans = ans+'H';
    		//else if(ss[i] == 'S' | ss[i] == 'G' | ss[i] == 'I' | ss[i] == 'T')
    		//  ans = ans+'A';//TODO - add this !
    			else
    				ans = ans+'C';

    		//System.out.println(resNum[i]+ "   " + aa[i] + "   " + ss[i] );
    	}
    	return ans;
    }

    public void readDSSP()  {
    	int i=0;
    	boolean cont = true;
    	String line;
    	try {
    		FileReader fr = new FileReader(fileName);
    		BufferedReader fdssp = new BufferedReader(fr);
    		line = fdssp.readLine();
    		while ((line != null) && (cont))
    		{
    			if (line.length() > 25)
    				if (line.substring(0,25).compareTo("  #  RESIDUE AA STRUCTURE") == 0)
    				{
    					cont = false;
    				}
    				else{
    					String word = line.substring(50,65);
    					if(word.equalsIgnoreCase(" Parallel BRIDG"))
    						parallelPrecentage  = Double.parseDouble(line.substring(6,10));
    					else if(word.equalsIgnoreCase("iParallel BRIDG"))
    						antyparallelPrecentage = Double.parseDouble(line.substring(6,10));
    					else if(word.equalsIgnoreCase("E O(I)-->H-N(J)") )
    						hbPrecentage = Double.parseDouble(line.substring(6,10));
    				}
    			line = fdssp.readLine();
    		} 
    		while (line != null)
    		{
    			i++;
    			line = fdssp.readLine();
    			if((line != null) && (line.indexOf("!")) > -1) i--;
    		} 		
    		fdssp.close();
    		fr = new FileReader(fileName);
    		fdssp = new BufferedReader(fr);
    		resNum = new int[i];
    		ss = new char[i];
    		aa = new char[i];
    		chainID = new char[i];
    		solvACC = new double[i];
    		for (int cc=0 ; cc<resNum.length ; cc++)
    			resNum[cc] = -999;
    		i=0;	
    		cont = true;
    		line = fdssp.readLine();
    		while ((line != null) && (cont))
    		{
    			if (line.length() > 25)
    				if (line.substring(0,25).compareTo("  #  RESIDUE AA STRUCTURE") == 0)
    				{
    					cont = false;
    				}
    			line = fdssp.readLine();
    		} 
    		while (line != null) {
    			if (line.indexOf("!") > -1) {
    				i--;
    			}
    			else {
    				resNum[i] = (Integer.valueOf(line.substring(5, 10).trim()));
    				aa[i] = line.charAt(13);
    				chainID[i] = line.charAt(11);
    				if (aa[i]>='a') // there is a funny way of marking Cystines in DSSP by lower case letters depending of their S-S indexing
    					aa[i] = 'C';
    				ss[i] = line.charAt(16);
    				if (ss[i] == ' ')
    					ss[i] = 'C';
    				solvACC[i] = (Double.valueOf(line.substring(35, 38).trim()));
    			}
    			i++;
    			line = fdssp.readLine();
    		} 		
    		fdssp.close();					
    	} // of the try
    	catch (Exception e) {
    		throw new RuntimeException(e);
    	} // of the catch
    }

    public SequenceList sequenceList() {
    	SequenceList seqList = new SequenceList();
    	Sequence aaSeq = new ResidueSequence(getAA(),"aa sequence of "+fileName);
    	seqList.add(aaSeq);
    	Sequence ssSeq = new SecondaryStructureSequence(getSSInOneLineConventional() ,"ss sequence of "+fileName);
    	seqList.add(ssSeq);
    	Sequence accSeq = new AccesibilitySequence(getRelativeACC() ,"accesibility sequence of "+fileName);                seqList.add(accSeq);
    	seqList.add(accSeq);
    	return seqList;
    }
}


/*
 * The dssp line:
          1         2         3         4         5         6         7         8         9         0
01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    3   13 A T  E >>  -A  399   0A  63     -2,-0.3     4,-1.1   396,-0.2     3,-0.7  -0.818  45.3 -28.2-122.9 157.9  -14.8   39.4   -6.4
 6241  520 N I              0   0   54   -743,-3.0  -741,-1.9    -2,-0.5    -2,-0.0  -0.899 360.0 360.0-115.3 360.0   33.5   14.7   14.2
*/
