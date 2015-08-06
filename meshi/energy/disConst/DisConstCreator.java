package meshi.energy.disConst;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.MeshiProgram;
import meshi.util.filters.Filter;


public class DisConstCreator extends EnergyCreator  implements KeyWords {

	private double numConstPerAtom = 2;
	private double maxDisForConst = 6; //Ang
	Filter filter = null;
	public DisConstCreator() {
		super(1.0);
	}

	public DisConstCreator(double weight, double maxDisForConst, double numConstPerAtom,  Filter filter) {
		super(weight);
		this.filter = filter;
		this.numConstPerAtom = numConstPerAtom;
		this.maxDisForConst = maxDisForConst;
	}


	public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
			CommandList commands) {
		return new  DisConstEnergy(createConstList(protein.atoms()), weight()/numConstPerAtom);
	}
	
	private DistanceList createConstList(AtomList atoms) {
		DistanceList list = new DistanceList();	   
		AtomList tmpList = null;
		double dis=0;
		int rand=0;
		for (int loop1=0 ; loop1<atoms.size() ; loop1++)
		if (filter.accept(atoms.atomAt(loop1))){
			tmpList = new AtomList();
			for (int loop2=loop1+1 ; loop2<atoms.size() ; loop2++)
			if (filter.accept(atoms.atomAt(loop2)) && 
						(Math.abs(atoms.atomAt(loop1).residueNumber()-atoms.atomAt(loop2).residueNumber())>4) ) {
				dis = atoms.atomAt(loop1).distanceFrom(atoms.atomAt(loop2)); 
				if (dis<maxDisForConst) {
					tmpList.add(atoms.atomAt(loop2));
				}
			}
			if (tmpList.size()>0) {
				for (int constCounter = 0; constCounter < numConstPerAtom ; constCounter++) {
					rand = (int)(MeshiProgram.randomNumberGenerator().nextDouble()*tmpList.size());
					list.add(new Distance(atoms.atomAt(loop1),tmpList.atomAt(rand)));
					//System.out.println("ADDed: " + rand + " " + tmpList.size() + " " + 
					//		atoms.atomAt(loop1).distanceFrom(tmpList.atomAt(rand)) + "\n" +
					//		atoms.atomAt(loop1) + "\n" + tmpList.atomAt(rand) +"\n");
				}
			}
			else {
				//System.out.println("List is empty for this atom:\n" + atoms.atomAt(loop1) + "\n");
			}
		}
		System.out.println(list.size() + " constraints created");
		return list;
	}
}
