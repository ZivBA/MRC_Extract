package meshi.util;
import java.util.Iterator;

/**
 * Key is a semantic sugar for a String that is used as a key to the ComandList class. 
 * It is intended to prevent confusion between the string we need (say file name) and 
 * the string used as a key to fish it from the commands list.
 **/
public class Key {
    public final String key;
    private static MeshiList keys = new MeshiList();
    public Key (String key) {
	this.key = key;
	for (Iterator iter = keys.iterator(); iter.hasNext();)
	    if( ((Key) iter.next()).key.equals(key))
		throw new RuntimeException("Key "+key+" is defined more than once.");
	keys.add(this);
    }
    
}

