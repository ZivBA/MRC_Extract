package meshi.energy;
import java.util.Arrays;
import java.util.Iterator;

import meshi.util.MeshiList;
import meshi.util.file.MeshiLineReader;
import meshi.util.filters.Filter;

/**
 * A list of parameters for an energy term.
 * Energy terms (extensions of AbstractEnergy) typically need a large number of parameters. These 
 * parameters are provided by the user in a text file. ParametersList objects read, parse and stores
 * the contennts of these files. They provide the parameters to the energy term with the 
 * getParameters(Parameters key) method. In general, each energy term requires a specific 
 * ParametersList class. See for example meshi.energy.bond.BondParametersList.
 **/
public abstract class ParametersList extends MeshiList {
    /**
     * True if the list is sortable.
     **/ 
    boolean sortable;

    /**
     * True if the list is sorted.
     **/
    boolean sorted = false;

    /**
     * Construct an empty list.
     **/
    public ParametersList() {
	super(new IsParameters());
    }

    /**
     * Construct a ParametersList object from parameters file. It is assumed that the file is arranged 
     * in raws where each raw includes the parameters for a specific interaction. the # sign may appear 
     * anywhere in a line and indicates the beginning of a comment. If the parameters are sortable the list 
     * is sorted. See for example the parameters for 
     * bond-energy in meshi/parameters/meshiPotential/bondParameters.dat .
     **/
    public ParametersList(String parametersFileName, boolean sortable) {
	super(new IsParameters());
	this.sortable = sortable;
	if (parametersFileName != null) {
	    System.out.println("Loading "+this+" parameters from "+parametersFileName);
	    try {
		MeshiLineReader lines = new MeshiLineReader(parametersFileName);
		String line;
		while ((line = lines.readLine("#")) != null) {
		    add(createParameters(line));
		}
		if (sortable) sort();
	    }
	    catch (RuntimeException e) {
		System.out.println("A problem while reading parameters file "+
				   parametersFileName);
		throw e;
	    }
	}
   }


    /** 
     * Construct a ParametersList object from multiple files.     
     * This constructor is useful when the parameters are distributed over many files,
     * where each file holds a single parameter. An example for this case are the TwoTorsions
     * energies.
     **/
    public ParametersList(String[] parametersFileName, boolean sortable) {
	super(new IsParameters());
	this.sortable = sortable;
	if (parametersFileName != null) {
    if (parametersFileName.length > 0) { 
	    System.out.println("Loading "+this+" parameters");
	    for (int cc=0; cc<parametersFileName.length ; cc++) {
	    	 add(createParameters(parametersFileName[cc]));
    		 if (sortable) sort();
 	    }
	}
	}
    }

      /**
      * Construct an empty list with Filter filter.
      **/
     public ParametersList(Filter filter) {
       super(filter);
     }

    /**
     * Adds an element (must be an instance of Parameters) to the list.
     **/   
    public boolean add(Object element) {
	if (sorted) throw new RuntimeException("Cannot add to ParametersList after sorting");
	return super.add(element);
    }
    
    /**
     * Sort the list. Note that if the List is not sortable this method simply does nothing.
     **/
    public void sort() {
	if (sortable) {
	    trim();
	    Arrays.sort(internalArray);
	    sorted = true;
	}
    }

    public Iterator iterator() {
	if (sortable & (! sorted)) 
	    throw new RuntimeException("Sortable ParametersList not sorted yet.\n"+
				       "Cannot provide an iterator");
	return super.iterator();
    }

    /** 
     * <pre>
     * Fetches a parameter from the list. A sorted list is assumed and a binary search is performed.
     * If the list is not sortable this method needs to be overrun. Example (from BondParametersList):
     *  
     * 	Parameters key = new BondParameters(pair.largeType(), pair.smallType());
     *  return getParameters(key);
     *</pre>
     **/ 
    public Parameters getParameters(Parameters key) {
	if (sortable) {
	    int index = Arrays.binarySearch(internalArray, key);
	    if (index <0) return null;
	    return (Parameters) elementAt(index);
	}
	throw new RuntimeException("The generic getParameters(Parameters key, \n"+
				   "                          ParametersList parametersList)\n"+
				   "method of ParametersList assume sorable parameters and \n"+
				   "uses binary search. Apparently this approach is not \n"+
				   "applicable here.");
    }

    /**
     * Energy term specific method to fetch parameters for the interactions between a st of atoms.
     **/
    public abstract Parameters parameters(Object Obj);


    /**
     * Energy term specific method to create a Parameters object from a line of the parameters file.
     **/
    public abstract Parameters createParameters(String line); 

    /**
     * Get an element from the list.
     **/
    public Parameters parametersAt(int index) {
	return (Parameters) elementAt(index);
    }

    //---------------------------------------------------------------------------------------------
    private static class IsParameters implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof Parameters);
	}
    }
}

    
    
    
