package meshi.energy;
/**
 * The parameters needed for the evaluation of an energy elements. 
 * There is not much in this class as the requirements of each energy 
 * term are very different from the requirements of the others.
 * See BondParameters for example.
 **/
public class Parameters {
    /** 
     * Converts a string to an int.
     **/
    public static int toInt(String s) {
	return Integer.valueOf(s.trim()).intValue();
    }
    /** 
     * Converts a string to a double.
     **/
    public static double toDouble(String s) {
	return Double.valueOf(s.trim()).doubleValue();
    }
}
    
