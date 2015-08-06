package meshi.energy.hydrogenBondsPlane;

import meshi.util.MeshiAttribute;

public class CN_AtomAttribute implements MeshiAttribute{
    public final boolean isC;
    public final boolean isN;
    public static final int key = CN_ATTRIBUTE;
    public CN_AtomAttribute(boolean isC, boolean isN) {
	if ((isC & isN) | (!((isC | isN))))
	    throw new RuntimeException("Weird parameters "+isC+" "+isN);
	this.isN = isN;
	this.isC = isC;
    }
    public final int key(){return key;}
    public String toString() {
	    	return "CN_AtomAttribute "+key+" "+isC+" "+isN;
    }
}
