package utils.UtilExceptions;

/**
 * Created by zivben on 03/10/16.
 */
public class ResidueLengthMismatch extends Throwable {
	
	public ResidueLengthMismatch(String name, String name1) {
		super("Tried to compare res: "+name+" with res: "+name1+" but they are of different lengths");
	}
}
