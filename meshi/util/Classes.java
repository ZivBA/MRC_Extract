package meshi.util;
import meshi.energy.LennardJones.LennardJonesEnergyElement;
import meshi.energy.LennardJones.LennardJonesParameters;
import meshi.energy.bond.BondEnergyElement;
import meshi.energy.bond.BondParameters;

public interface Classes {
    Class BOND_PARAMETERS_CLASS = (new BondParameters()).getClass();    
    Class BOND_ELEMENT_CLASS    = (new BondEnergyElement()).getClass();
    Class LENNARD_JONES_PARAMETERS_CLASS = (new LennardJonesParameters()).getClass();    
    Class LENNARD_JONES_ELEMENT_CLASS    = (new LennardJonesEnergyElement()).getClass();
}
