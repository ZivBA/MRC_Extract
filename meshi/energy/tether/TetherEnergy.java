package meshi.energy.tether;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.MeshiList;


public class TetherEnergy  extends SimpleEnergyTerm{
    /**
     * The constructor associates any bond with its parameters.
     **/
    private int atomCounter=0;

    public TetherEnergy() {}

    public TetherEnergy(AtomList atomList, AtomList takePegFrom,  double weight) {
	
	elementsList = new MeshiList();
	this.weight = weight;

	for (atomCounter=0 ; atomCounter<atomList.size(); atomCounter++) {
		Atom atom = atomList.atomAt(atomCounter);
		if (takePegFrom!=null){
			if (atom.reliability()>0.0) {
				Atom pegAtom = takePegFrom.findAtomInList(atom.name(), atom.residueNumber());
				if (pegAtom==null)
					System.out.println("Not taking constraints to: " + atom.name() + atom.residueNumber() + "because they are not in the file you specified");
				else {
					EnergyElement newElement = createTetherElement(atom,pegAtom);
					if (!newElement.frozen())
						elementsList.add(newElement);
				}
			}
		}
		else {
			if (atom.reliability()>0.0) {
				EnergyElement newElement = createTetherElement(atom,null);
				if (!newElement.frozen())
					elementsList.add(newElement);
			}			
		}
	}
	comment = "Tether";
    }


    public EnergyElement createTetherElement(Atom baseElement , Atom peg) {
    	return new TetherEnergyElement(weight, baseElement , peg);
    }

     public void update(){}
     public void update(int i ){}

     public EnergyElement createElement(Object baseElement,Parameters p){return null;}
     
     /**
      * The pegs will be set to the current position of the atoms.
      */
     public void updatePegsToCurrentPosition() {
    	 for (int cc=0; cc<elementsList().size() ; cc++) {
    		 TetherEnergyElement element = (TetherEnergyElement) elementsList().elementAt(cc);
    		 element.setNewPegCoordinates(element.atom.x(), element.atom.y(), element.atom.z());
    	 }
     }

}


