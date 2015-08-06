package meshi.PDB;
import java.util.Formatter;

import meshi.util.MeshiException;
/**

from http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html

   COLUMNS        DATA TYPE       FIELD         DEFINITION
---------------------------------------------------------------------------------
 1 -  6        Record name     "ATOM  "

 7 - 11        Integer         serial        Atom serial number.

13 - 16        Atom            name          Atom name.

17             Character       altLoc        Alternate location indicator.

18 - 20        Residue name    resName       Residue name.

22             Character       chainID       Chain identifier.

23 - 26        Integer         resSeq        Residue sequence number.

27             AChar           iCode         Code for insertion of residues.

31 - 38        Real(8.3)       x             Orthogonal coordinates for X in
                                             Angstroms.

39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in
                                             Angstroms.

47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in
                                             Angstroms.

55 - 60        Real(6.2)       occupancy     Occupancy.

61 - 66        Real(6.2)       tempFactor    Temperature factor.

73 - 76        LString(4)      segID         Segment identifier, left-justified.

77 - 78        LString(2)      element       Element symbol, right-justified.

79 - 80        LString(2)      charge        Charge on the atom.
 */
public class PdbLine  {
    private String line;
    public PdbLine(String line) {
	this.line = line;
    }
    public  PdbLine (int number, String name, String residueName,
                     String chain, int residueNumber,
                     double x,double y, double z, double occupancy,
                     double temperatureFactor) {
        StringBuffer res = new StringBuffer();
        Formatter f = new Formatter(res);

        // occupancy is between 0.00 and 1.00
        if(occupancy < 0)
            occupancy = 0;
        else if(occupancy > 1)
            occupancy = 1;
        // temperature factor's threshold are -99 and 999
	if(temperatureFactor<=-99)
            temperatureFactor = 99.0;
	else if(temperatureFactor>=999.0)
            temperatureFactor = 999.0;

        switch(name.length()){
        case 1:name = " "+name+"  ";break;
        case 2:name = " "+name+" ";break;
        case 3:name = " "+name;break;
        case 4:name = name;break;
        }
        
        f.format("%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s",
                 "ATOM",
                 number,
                 name,
                 " ",
                 residueName,
                 chain,
                 residueNumber,
                 " ",
                 x,y,z,
                 occupancy,temperatureFactor,
                 " "," "," ");
        f.flush();
        f.close();
        line = res.toString();
    }
    
    public boolean isAnAtom() {
	return line.startsWith("ATOM");
    }
    public boolean isAHeteroAtom() {
	return line.startsWith("HETATM"); 
    }
    public boolean isAnAtomOrHeteroAtom() {
	return (isAnAtom() || isAHeteroAtom());
    }
    public void needsToBeAnAtom() {
	if (! isAnAtomOrHeteroAtom() )
	    throw new MeshiException("PdbLine error:\n"+
                                     "needs to be an atom:\n"+
                                     this+"\n");
    }
    public boolean isAComment() {return (! isAnAtomOrHeteroAtom());}
    public boolean isSEQRES() {return line.startsWith("SEQRES");}
    public String toString() {return line;}
    public double x() {
	try {
	    needsToBeAnAtom();
	    return Double.valueOf(line.substring(30,38).trim()).doubleValue();
	}
	catch (Exception ex) { 
	    throw new RuntimeException("PdbLine Error - Badly formatted PdbLine:\n"+line+"\n"+ex);
	}
    }
    public double y() {
	try {
	    needsToBeAnAtom();
	    return Double.valueOf(line.substring(38,46).trim()).doubleValue();
	}
	catch (Exception ex) { 
	    throw new RuntimeException("PdbLine Error - Badly formatted PdbLine:\n"+line+"\n"+ex);
	}
     }
    public double z() {
	try {
	    needsToBeAnAtom();
	    return Double.valueOf(line.substring(46,54).trim()).doubleValue();
	}
	catch (Exception ex) { 
	    throw new RuntimeException("PdbLine Error - Badly formatted PdbLine:\n"+line+"\n"+ex);
	}
    }    
    public String chain(){
	needsToBeAnAtom();
	String temp = line.substring(21,22);
	//	if (temp.equals(" ")) temp = "0";  commented out 28.8.2002 by chen
	return temp;
    }
    public String residueName() {
	needsToBeAnAtom();
	return line.substring(17,20).trim();
    }
    public String name() {
	return line.substring(12,16).trim();
    }
    public Integer residueNumber() {
	needsToBeAnAtom();
	return Integer.decode(line.substring(22,26).trim());
    }
    public int number() {
	needsToBeAnAtom();
	return Integer.valueOf(line.substring(6,11).trim()).intValue();
    }
    /*
     *Check if this is a MODEL line.
     */
    public boolean isAModel() {
	return line.startsWith("MODEL");
    }    
    /*
     *If this is a MODEL line, the function returns the model number, otherwise it returns -1.
     */
    public int getModel() {
	if (isAModel()) 
            return Integer.valueOf(line.substring(6).trim()).intValue();
	else
            return -1;
    }
    public double temperatureFactor() {
	needsToBeAnAtom();
	if (line.length() > 65)
            try{return Double.valueOf(line.substring(60,66).trim()).doubleValue();}
            catch(NumberFormatException nfe){return 0.0;}
	else
            return 0.0;
    }     
    public double occupancyNumber() {
	needsToBeAnAtom();
	if (line.length() > 59)
            try{return Double.valueOf(line.substring(54,60).trim()).doubleValue();}
            catch(NumberFormatException nfe){return 0.0;}
	else
            return 0.0;
    }           
    public String alternateLocation() {
	needsToBeAnAtom();
	return line.substring(16,17).trim();
    }   
}
