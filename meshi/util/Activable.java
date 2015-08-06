package meshi.util;

/**
 * In order to allow inactivation/activation of atoms (essentially the same as freezing
 * that do not kill the distance matrix), this interface should be implemented by any Activable object
 * that is benefiting from the inactivation procedure. Any such element must report to the atom
 * it depends on about its existence by the use of the Atom.addActivableToList(this) method. 
 * When an atom activation status is changed, it will invoke the 'updateActivity' method
 * in that element.
 * 
 * The method 'addObjectToActivableListInAtoms' should arrange for the addition of the object to 
 * the relevant atoms. It should invoked from the constructor, once the atoms involved in the object are known.
 */

public interface Activable {
	public void updateActivity();	
	public void addObjectToActivableListInAtoms();
	public boolean active();
}
