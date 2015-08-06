/*
 * Created on 16/11/2004
 * Window - Preferences - Java - Code Style - Code Templates
 */
package meshi.energy.hydrogenBond;

import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.energy.LennardJones.LennardJonesParametersList;
import meshi.geometry.Distance;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;

/**
 * @author amilev
 *
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class HydrogenBondsParametersList extends LennardJonesParametersList {
	AtomList nitrogens;
	
	public HydrogenBondsParametersList(String parametersFileName,AtomList nitrogens){
		super(parametersFileName);
		this.nitrogens = nitrogens;
	}

    public HydrogenBondsParametersList(String parametersFileName){
		super(parametersFileName);
    }
	
	/* (non-Javadoc)
	 * @see meshi.energy.ParametersList#getParameters(java.lang.Object)
	 */
	public Parameters parameters(Object obj) {
		Distance distance = (Distance) obj;
		int largeType = distance.largeType();
		int smallType = distance.smallType();
		try {
			return (HydrogenBondsParameters) elementAt(largeType*(largeType+1)/2+smallType);
	}	
		catch (Exception e) {
			throw new RuntimeException(largeType+"-"+Atom.type(largeType)+" "+ 
					smallType+"-"+Atom.type(smallType)+"\n"+e); 
		}	
	}

	/* (non-Javadoc)
	 * @see meshi.energy.ParametersList#createParameters(java.lang.String)
	 */
	public Parameters createParameters(String line) {
			return new HydrogenBondsParameters(new StringTokenizer(line));
	}

	/**
	 * @return nitrogens
	 */
	public AtomList nitrogens() {
		return nitrogens;
	}

}
