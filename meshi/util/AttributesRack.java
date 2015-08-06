package meshi.util;

public class AttributesRack{
    public static final int MAX_ATTRIBUTES = 10;
    private final Object[] internalArray = new Object[MAX_ATTRIBUTES];
    public void addAttribute(MeshiAttribute attribute){
	int key = attribute.key();
	if (key > MAX_ATTRIBUTES) throw new RuntimeException("Key = "+key+" is larger than "+
							     "MAX_ATTRIBUTES "+MAX_ATTRIBUTES);
	internalArray[key] = attribute;
    }
    public MeshiAttribute getAttribute(int key) {
	return (MeshiAttribute) internalArray[key];
    }
}
 
