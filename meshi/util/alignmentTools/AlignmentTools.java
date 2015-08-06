package meshi.util.alignmentTools;
import java.util.Iterator;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.KeyWords;


public class AlignmentTools implements Residues, AtomTypes , KeyWords { 

	public static void thread(Protein prot, int begins, int ends, int move) {
		AtomList a1,a2;
		Atom at1,at2;
		for (int res=begins ; res<=ends ; res++) {
			a1 = prot.residue(res).atoms().filter(new AtomList.BackboneFilter());
			a2 = prot.residue(res+move).atoms().filter(new AtomList.BackboneFilter());
			for (int c=0 ; c<a1.size() ; c++) {
				at1 =a1.atomAt(c);
				at2 = getAtom(a2,at1.name());
				if (at2!=null) 
					at1.setXYZ(at2.x(),at2.y(),at2.z());
			}
		} 
	}
	
	public static void outOfTheWay(Protein prot, int resnum, double shift)	{
		double x,y,z,cmx=0, cmy=0, cmz=0; // center of mass x, y and z
    	Atom atom;
		Iterator iter = prot.atoms().iterator(); 
		while((atom = (Atom) iter.next()) != null) { 
	    	cmx += atom.x();
		    cmy += atom.y();
		    cmz += atom.z();
		}
		cmx /= prot.atoms().size();
		cmy /= prot.atoms().size();
		cmz /= prot.atoms().size();
		
		AtomList a1 = prot.residue(resnum).atoms();
		for (int c=0 ; c<a1.size() ; c++) {
			x=a1.atomAt(c).x();
			y=a1.atomAt(c).y();
			z=a1.atomAt(c).z();
			double norm = Math.sqrt((x-cmx)*(x-cmx) + 
									(y-cmy)*(y-cmy) +
									(z-cmz)*(z-cmz));
			a1.atomAt(c).setX(cmx + (norm+shift)/norm*(x-cmx));
			a1.atomAt(c).setY(cmy + (norm+shift)/norm*(y-cmy));
			a1.atomAt(c).setZ(cmz + (norm+shift)/norm*(z-cmz));
		}
	}	
	
	public static Atom getAtom(AtomList al , String name) {
		for (int c=0 ; c<al.size() ; c++) 
			if (al.atomAt(c).name().equals(name))
				return al.atomAt(c);
		return null;
	}

	public static Atom getAtom(AtomList al , int resNum, String name) {
		for (int c=0 ; c<al.size() ; c++) 
			if ((al.atomAt(c).residueNumber()==resNum) && al.atomAt(c).name().equals(name))
				return al.atomAt(c);
		return null;
	}
}




