package utils.UtilExceptions;

import utils.molecularElements.SimpleProtein;

/**
 * Created by zivben on 08/09/15.
 */
public class MissingChainID extends Throwable {
	public MissingChainID(char chainID, SimpleProtein tempProt) {
		super("Tried to get chain " + chainID + "\nAs requested by SCWRL protein: " + tempProt.getSource()
				.getAbsolutePath());
	}
}
