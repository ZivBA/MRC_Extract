package meshi.geometry;
import java.util.Iterator;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.AtomPairList;
import meshi.molecularElements.Protein;
import meshi.util.MeshiIterator;
import meshi.util.MeshiList;
import meshi.util.Updateable;
import meshi.util.UpdateableException;
import meshi.util.filters.Filter;
/**
 **/
public class TorsionList extends MeshiList implements Updateable {
     private int numberOfUpdates = 0;
   /**
     * An empty Torsion list
     **/
    public TorsionList() {
        this(new IsTorsion());
    }
    /**
     * An empty TorsionList
     **/
     protected TorsionList(Filter filter) {
        super(filter);
    }




    /**
     * An Torsion list based on a angle list
     **/
    public TorsionList(AngleList angles, DistanceMatrix distanceMatrix) {
	this();
	MeshiIterator angleIter = angles.meshiIterator();
	Torsion tor;
	Angle angle1, angle2;	
	while((angle1 = (Angle) angleIter.next()) != null) {
	    while((angle2 = (Angle) angleIter.nestedNextTo()) != null) 
		if ((angle1.sharedAtomPair(angle2) != null) &&      // There is a shared pair and the angles are not the same but in opposite directions
		    !((angle1.atom1==angle2.atom3)&&(angle1.atom2==angle2.atom2)&&
		      (angle1.atom3==angle2.atom1))) {
		    tor = getTorsion(angle1, angle2, distanceMatrix);
		    if (isNamed(tor))
			add(tor);
		    tor =  getTorsion(angle2, angle1, distanceMatrix);
		    if (isNamed(tor))
			add(tor);
		}
	}
    }

      public Torsion getTorsion(Angle angle1, Angle angle2, DistanceMatrix distanceMatrix) {
	     return new Torsion(angle1, angle2, distanceMatrix);
      } 

    public void update(int numberOfUpdates) throws UpdateableException {
	if (numberOfUpdates == this.numberOfUpdates+1) {
	    int size = size();
	    for (int i = 0; i < size; i++) {
		torsionAt(i).update(numberOfUpdates);
	    }
	    this.numberOfUpdates++;
	}
	else if (numberOfUpdates != this.numberOfUpdates) 
	    throw new RuntimeException("Something weird with TorsionList.update(int numberOfUpdates)\n"+
				       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
    }


    public Torsion torsionAt(int i) { return (Torsion) elementAt(i);}
    
    public AtomList atomList() {
	AtomList list = new AtomList();
	Iterator iter = iterator();
	Torsion torsion;
	Atom atom1, atom2, atom3, atom4;
	while ((torsion = (Torsion) iter.next()) != null) {
	    atom1 = torsion.atom1;
	    atom2 = torsion.atom2;
	    atom3 = torsion.atom3;
	    atom4 = torsion.atom4;
	    if (! list.contains(atom1)) {
		list.add(atom1);
	    }
	    if (! list.contains(atom2)) {
		list.add(atom2);
	    }
	    if (! list.contains(atom3)) {
		list.add(atom3);
	    }
	    if (! list.contains(atom4)) {
		list.add(atom4);
	    }
	}
	return list;
    }

    public boolean equivalentExists(Torsion findMe) {
	Iterator torsions = iterator();
	Torsion torsion;
	while ((torsion = (Torsion) torsions.next()) != null)
	    if (torsion.equivalent(findMe)) return true;
	return false;
    }
	
    public TorsionList filterEquivalents() {
	TorsionList out = new TorsionList();
	Iterator torsions = iterator();
	Torsion torsion;
	while ((torsion = (Torsion) torsions.next()) != null) {
	    if (isPhi(torsion)) out.add(torsion);
	    else if (isPsi(torsion)) out.add(torsion);
	    else if (isOmega(torsion)) out.add(torsion);
	    else if (isChi1(torsion)) out.add(torsion);
	    else if (isChi2(torsion)) out.add(torsion);
	    else if (isChi3(torsion)) out.add(torsion);
	    else if (isChi4(torsion)) out.add(torsion);
	    else if (isNimp(torsion)) out.add(torsion);
	    else if (isCimp(torsion)) out.add(torsion);
	}
	torsions = iterator();
	while ((torsion = (Torsion) torsions.next()) != null) 
	    if (! out.equivalentExists(torsion)) out.add(torsion);
	return out;
    }
	
    /**
     * Returns a sub-list that is accepted by the parameter
     **/
    public TorsionList chi1Filter() {
	Iterator torsions = iterator();
	Torsion torsion;
	TorsionList out = new TorsionList();
	while ((torsion = (Torsion) torsions.next()) != null) 
	    if (isChi1(torsion)) out.add(torsion);
	return out;
    }

    /**
     * Returns a sub-list that is has a known name
     **/    
    public TorsionList namedFilter() {
	Iterator torsions = iterator();
	Torsion torsion;
	TorsionList out = new TorsionList();
	while ((torsion = (Torsion) torsions.next()) != null) 
	    if (isNamed(torsion)) out.add(torsion);
	return out;
    }

    public static boolean isNameInIUPAC(Torsion torsion) {
       	if (isPhi(torsion) ||
       			isPsi(torsion) ||
       			isOmega(torsion) ||
       			isChi1(torsion) ||
       			isChi2(torsion) ||
       			isChi3(torsion) ||
       			isChi4(torsion))       			
    	   return true;
    	else
    	   return false;
    }

    
    public static boolean isNamed(Torsion torsion) {
       	if (torsion.getTorsionName().compareTo("") != 0)
    	   return true;
    	else
    	   return false;
    }

    public static boolean isPhi(Torsion torsion) {
       	if (torsion.getTorsionName().compareTo("PHI") == 0)
    	   return true;
    	else
    	   return false;
    }
    
    public static boolean isPsi(Torsion torsion) {
       	if (torsion.getTorsionName().compareTo("PSI") == 0)
    	   return true;
    	else
    	   return false;
    }
    
    public static boolean isOmega(Torsion torsion) {
       	if (torsion.getTorsionName().compareTo("OMG") == 0)
    	   return true;
    	else
    	   return false;
    }
    
    public static boolean isCimp(Torsion torsion) {
	if (torsion.atom1.name.equals("CA") &
	    torsion.atom2.name.equals("N") &
	    torsion.atom3.name.equals("C") &
	    torsion.atom4.name.equals("O")) return true;
	return false;
    }
    
    public static boolean isNimp(Torsion torsion) {
	if (torsion.atom1.name.equals("H") &
	    torsion.atom2.name.equals("C") &
	    torsion.atom3.name.equals("N") &
	    torsion.atom4.name.equals("CA")) return true;
	return false;
    }

    public static boolean isChi1(Torsion torsion) {
    	if (torsion.getTorsionName().compareTo("CHI1") == 0)
    	   return true;
    	else
    	   return false;
    }
    public static boolean isChi2(Torsion torsion) {
    	if (torsion.getTorsionName().compareTo("CHI2") == 0)
    	   return true;
    	else
    	   return false;
    }
    public static boolean isChi3(Torsion torsion) {
    	if (torsion.getTorsionName().compareTo("CHI3") == 0)
    	   return true;
    	else
    	   return false;
    }
    public static boolean isChi4(Torsion torsion) {
    	if (torsion.getTorsionName().compareTo("CHI4") == 0)
    	   return true;
    	else
    	   return false;
    }

    public static boolean isOOP(Torsion torsion) {
    	if (torsion.getTorsionName().compareTo("OOP") == 0)
    	   return true;
    	else
    	   return false;
    }

    
    public static class IsTorsion implements Filter {
	public boolean accept(Object obj) {
		return (obj instanceof Torsion); 
        }
    }

    public static class FilterOOP implements Filter {
	public boolean accept(Object obj) {
		Torsion tor = (Torsion) obj;
		return isOOP(tor); 
        }
    }
    
    public static class FilterPhi implements Filter {
	public boolean accept(Object obj) {
		Torsion tor = (Torsion) obj;
		return isPhi(tor); 
        }
    }

    public static class FilterPsi implements Filter {
	public boolean accept(Object obj) {
		Torsion tor = (Torsion) obj;
		return isPsi(tor); 
        }
    }

    public static class FilterChi1 implements Filter {
	public boolean accept(Object obj) {
		Torsion tor = (Torsion) obj;
		return isChi1(tor); 
        }
    }

    public static class FilterChi2 implements Filter {
        public boolean accept(Object obj) {
                Torsion tor = (Torsion) obj;
                return isChi2(tor);
        }
    }

    public static class FilterChi3 implements Filter {
        public boolean accept(Object obj) {
                Torsion tor = (Torsion) obj;
                return isChi3(tor);
        }
    }

    public static class FilterChi4 implements Filter {
        public boolean accept(Object obj) {
                Torsion tor = (Torsion) obj;
                return isChi4(tor);
        }
    }

    public static class FilterSideChain implements Filter {
        public boolean accept(Object obj) {
                Torsion tor = (Torsion) obj;
                return (isChi1(tor) || isChi2(tor) || isChi3(tor) || isChi4(tor));
        }
    }

    public static class FilterFamouseTorsions implements Filter {
        public boolean accept(Object obj) {
                Torsion tor = (Torsion) obj;
                return (isChi1(tor) || isChi2(tor) || isChi3(tor) || isChi4(tor) ||
                		isPhi(tor) || isPsi(tor));
        }
    }

    public boolean sortable() {return false;}

    public static TorsionList createTorsionList(Protein protein, DistanceMatrix distanceMatrix) {
	AtomPairList bondList = protein.bonds();
	bondList.renumber();
	AngleList angleList = new AngleList(bondList, distanceMatrix);
	return new TorsionList(angleList, distanceMatrix);
    }
    
    public static TorsionList createQuickAndDirtyTorsionList(Protein protein, DistanceMatrix distanceMatrix) {
	    AtomPairList bondList = protein.bonds();
	    bondList.renumber();
	    AngleList angleList = new AngleList(bondList, distanceMatrix);
	    return new QuickAndDirtyTorsionList(angleList, distanceMatrix);
    }
    
    public void freeze() {
	Iterator iter = iterator();
	Torsion torsion;
	while ((torsion = (Torsion) iter.next()) != null) 
	   torsion.freeze();
    }
    
    
}

