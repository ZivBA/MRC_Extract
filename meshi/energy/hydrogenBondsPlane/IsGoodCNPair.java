package meshi.energy.hydrogenBondsPlane;


import meshi.geometry.Distance;
import meshi.util.filters.Filter;

public class IsGoodCNPair implements Filter{
    private IsCN isCN = new IsCN();
    private IsNot13 isNot13 = new IsNot13();
    private IsGoodSS goodSS = new IsGoodSS();
    	
    public boolean accept(Object obj) {
        Distance dis = (Distance)obj;
        if (dis.bonded()) return false;
        return isCN.accept(obj) && isNot13.accept(obj) && goodSS.accept(obj);
    }
 //--------------------------- internal class IsNot13 -----------------------------------
    /* this filter Accept Distance between Atoms just if there are at list 3 atoms away on the sequense,
     *  (atomes that are less then 3 atoms away on the sequense never create Hydrogen bonds with each other)
     */
    static class IsNot13 implements Filter{
        public boolean accept(Object obj) {
            Distance dis = (Distance)obj;
        	if (!dis.atom1().chain().equalsIgnoreCase(dis.atom2().chain()))
        		return true;
            int residueNumber1 = dis.atom1().residueNumber();
            int residueNumber2 = dis.atom2().residueNumber();
            return (!(Math.abs(residueNumber1-residueNumber2) <= 2));
        }
    }

        
    //--------------------------- internal class GoodSS ---------------------------
    
    static class IsGoodSS implements Filter{
    	public boolean accept(Object obj) {
            Distance dis = (Distance)obj;  
    		String atom1SS = dis.atom1().residue().secondaryStructure();
    		String atom2SS = dis.atom2().residue().secondaryStructure();
    		if (atom1SS.equals("COIL") || atom2SS.equals("COIL"))
    			return false;
            if ((atom1SS.equals("HELIX") || atom1SS .equals("HELIX_OR_COIL")) && (atom2SS.equals("HELIX") | atom2SS .equals("HELIX_OR_COIL") ))
                            {
                                int residueNumber1 = dis.atom1().residueNumber();
                                int residueNumber2 = dis.atom2().residueNumber();
                                return (Math.abs(residueNumber1-residueNumber2) == 4);
                            }
//            return !(((atom1SS.equals("HELIX") || atom1SS .equals("HELIX_OR_COIL")) && (atom2SS.equals("SHEET") | atom2SS .equals("SHEET_OR_COIL") )) ||
//                    ((atom1SS.equals("SHEET") | atom1SS .equals("SHEET_OR_COIL")) && (atom2SS.equals("HELIX") | atom2SS .equals("HELIX_OR_COIL") )));
                        return (atom1SS.equals("SHEET") || atom1SS .equals("SHEET_OR_COIL")) && (atom2SS.equals("SHEET") | atom2SS .equals("SHEET_OR_COIL"));
          		 }
    }
}
