package meshi.energy.tether;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.Key;
import meshi.util.KeyWords;
import meshi.util.filters.Filter;


public class TetherCreator extends EnergyCreator  implements KeyWords {
	Filter filter = null;
	AtomList takePegFrom = null;

	public TetherCreator() {
		super( TETHER_ENERGY);
	}

	
	/**
	 * With this creator the tether's "peg" is taken from the model given to the createEnergyTerm by the TotalEnergy.
	 **/
	public TetherCreator(double weight, Filter filter) {
		super(weight);
		this.filter = filter;
	}

	public TetherCreator(Key key) {
		super(key);

	}


	/**
	 * If you invoke this method the tether's "peg" is taken from an atom in file: "fileName"
	 **/
	public void takePegFrom(String fileName) {
		if (fileName == null)
			takePegFrom = null;
		else
			takePegFrom = new AtomList(fileName);
	}
	
	
	public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
			CommandList commands) {
		AtomList atomList = protein.atoms();
		for (int j=0;j<atomList.size();j++) {
			if (filter != null) {
				if (filter.accept(atomList.atomAt(j))) 
					atomList.atomAt(j).setReliability(1);
				else 
					atomList.atomAt(j).setReliability(0);
			}
			else atomList.atomAt(j).setReliability(1);
		}

		return new  TetherEnergy(atomList, takePegFrom , weight());
	}
}
