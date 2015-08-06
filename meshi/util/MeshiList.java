package meshi.util;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;

import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;

/**
 * Type (class/interface) specific list with nested loop iterator.
 * Much of Structural biology software has to do with 
 * looping over large arrays with the following characteristics: <br>
 * <ol>
 *    <li> All objects are of the same type (e.g. All atoms of 
 *         the ...., All residues of....., etc. ). 
 *    <li> In such "uniform type" list <b>null</b> is not expected.
 *    <li> Nested loops are common.
 *    <li> The lists are generated once in some initialization phase of the program 
 *         and remain structurally unchanged during execution. Atoms of a protein for example
 *         are not added or removed during simulation (they may of course change their position
 *         velocity etc.).
 * </ol>
 * In MeshiList we tried to enforce and utilize this characteristics.
 * <p>
 * <ul>
 *    <li> The Iterator over MeshiList utilizes the fact that <b>null</b> can never be 
 *         an element of the list in the following way: next() returns <b>null</b> if and only if
 *         hasNext() == false.
 *    <li> MeshiList also offer another non-standard type of iterator 
 *         ({@link meshi.util.MeshiIterator MeshiIterator}) that makes it easier to write nested 
 *         loops over the list.
 * </ul> 
 **/
public class MeshiList {
    protected Object[] internalArray;
    protected int capacity;
    protected int size;
    protected int modCount;
    private static final int DEFAULT_CAPACITY = 500;
    public static final Filter isMeshiListFilter = new IsMeshiListFilter();
    /**
     * Describes the list.
     */
    private String comment;
    /**
     * Filters out unwanted objects
     */
    private Filter filter;

    public MeshiList(int capacity) {
    	this(new KolDichfin(), capacity);
    }
    public MeshiList() {
    	this(new KolDichfin());
    }
    /**
     * Creates an empty list.
     */
    public MeshiList(Filter filter) { 
	if (filter == null) throw new RuntimeException("filter is null");
	internalArray = new Object[DEFAULT_CAPACITY];
	capacity = DEFAULT_CAPACITY;
	size = 0;
	modCount = 0;
	this.filter = filter;
	
    }
    /**
     * Creates an empty list.with length = capacity
     */
      public MeshiList(Filter filter, int capacity) { 
	this.capacity = capacity;
        internalArray = new Object[capacity];       
	size = 0;
	modCount = 0;
	this.filter = filter;	
    }

    /**
     * Creates a list with the elements of the array parameter;
     **/ 
    public MeshiList(Filter filter, Object[] elements) { 
    	this(filter,elements,false);	
    }	

    /**
     * Creates a list with the elements of the array parameter and alow to ignor elements from 
     * differents types only is exclusiveAdd is true
     **/ 
    public MeshiList(Filter filter, Object[] elements, boolean exclusiveAdd){
    	this(filter);
    	int length = elements.length;
    	for (int i = 0; i < length; i++)
    	 	if(! exclusiveAdd)
    	 		add(elements[i]);
    	 	else
    	 	 	exclusiveAdd(elements[i]);
    }	

    
    public MeshiList(MeshiListIterator iterator) throws Exception{
        this(iterator.listFilter());
        Object obj;
        while ((obj=iterator.next())!=null) add(obj);
    }
   

    /** 
     * Adds only objects of the predefined type.
     **/      
    public boolean add(Object element) {
	if (! filter.accept(element)) {
	    throw new MeshiException("Cannot add "+element+" of class "+
			    element.getClass()+" to "+this);
	}
 	if (size < capacity) {
	    internalArray[size] = element;
	    size++;
	    modCount++;
	    return true;
	}
	else {
	    capacity *= 2;
	    Object[] newArray = new Object[capacity];
	    for (int i = 0; i < size; i++)
		newArray[i] = internalArray[i];
	    internalArray = newArray;
	    return add(element) ;
	}
    }
    
    public void add(MeshiList list, boolean breakFlag) {
	if (breakFlag) {
           Iterator iter = list.iterator();
	   Object obj;
	   while((obj = iter.next()) != null) add(obj);
	}
	else {
		if (! filter.accept(list)) {
			throw new MeshiException("Cannot add "+list+" of class "+
						 list.getClass()+" to "+this);
		}
		if (size < capacity) {
			internalArray[size] = list;
			size++;
			modCount++;
		}
		else {
			capacity *= 2;
			Object[] newArray = new Object[capacity];
			for (int i = 0; i < size; i++)
				newArray[i] = internalArray[i];
			internalArray = newArray;
			add(list, breakFlag);
		}
	}
    }
			

  
    /** 
     * Adds objects of the predefined type and ignor others.
     * @return true if element was added
     **/      
    public boolean exclusiveAdd(Object element) {
	if (! filter.accept(element)) {
	   return false;
	}
 	if (size < capacity) {
	    internalArray[size] = element;
	    size++;
	    modCount++;
	    return true;
	}
	else {
	    capacity *= 2;
	    Object[] newArray = new Object[capacity];
	    for (int i = 0; i < size; i++)
		newArray[i] = internalArray[i];
	    internalArray = newArray;
	    return add(element) ;
	}
    }
        
    /*
      * @param element
      * @return true
      * NOTE: this method doesn't note check that the element is approved to be added.
      */
    public boolean fastAdd(Object element) {
       if (size < capacity) {
           internalArray[size] = element;
           size++;
           modCount++;
           return true;
       }
       else {
           capacity *= 2;
           Object[] newArray = new Object[capacity];
           for (int i = 0; i < size; i++)
               newArray[i] = internalArray[i];
           internalArray = newArray;
           return fastAdd(element) ;
       }
     }

    public void add(MeshiList other) {
	    add(other,true);
    }

    public void set(int index, Object element) {
	if (! filter.accept(element)) {
	    throw new MeshiException("Cannot add "+element+" of class "+
			    element.getClass()+" to "+this);
	}
	if (index < capacity) {
	    internalArray[index] = element;
	    modCount++;
	}
	else {
	    throw new MeshiException("index out of bounds "+index);
	}	 
    }

    public Object remove(int index) {
	if (index >= size) 
	    throw new RuntimeException("index out of bounds "+index);
	Object out = internalArray[index];
	for (int i = index; i < size - 1; i++)
	    internalArray[i] = internalArray[i+1];
	size--;
	modCount++;
	return out;
    }
   

    public final Object elementAt(int index) {
	if (index >= size) 
	    throw new RuntimeException("index "+index+" out of bounds "+size+" "+internalArray.length);
	return internalArray[index];
    }

    public final Object fastElementAt(int index) {
	return internalArray[index];
    }
    
    public boolean contains(Object obj) {
	for(int i = 0; i < size; i++)
	    if (obj.equals(internalArray[i])) return true;
	return false;
    }


    /**
     * Returns a {@link meshi.util.MeshiIterator MeshiIterator} over the list elements
     **/
    public MeshiIterator meshiIterator() {
	return new MeshiIterator(this);
    }
    /**
     * Returns an Iterator over the list objects. 
     * In the Iterator interface the value returned by next() when hasNext() == false
     * is undefined. In this implementation next() returns <b>null</b> if and only if
     * hasNext() == false.
     **/
    public Iterator iterator() {
	Iterator out = new MeshiListIterator();// MeshiListIterator is an internal class.
	if (((MeshiListIterator) out).listFilter() == null) throw new RuntimeException("filter is null");
	return out; 
    }      
   
    public Iterator reversIterator() {
    	return (Iterator) new ReversIterator();
    } 
    // New methods not found in Vector
    /** 
     * Set the list comment.
     **/
    public void setComment(String s) {
	comment = s;
    }

    /**
     * Returns the list comment.
     **/
    public String comment() {return comment;}
    
    /**
     * Returns the internal filter of the list.
     */
    public Filter filter() {return filter;}
    

    public Object[] toArray() {
	Object[] out = new Object[size];
	for (int i = 0; i < size; i++) {
	    out[i] = internalArray[i];
	}
	return out;
    }

    /**
     * Returns an array of the list elements sorted by <b>comp</b>.
     **/
    public Object[] sortToArray(Comparator comp) {
	Object[] array = toArray();
	Arrays.sort(array,comp);
	return array;
    }

    /**
     * Print the list elements.
     **/
    public void print() {
	print(1);
    }
    /**
     * Print the list elements in raws. 
     * The number of elements in a raw is determined by the parameter
     **/
    public void print(int rawLength) {
	if (rawLength < 1) throw new RuntimeException("Raw Length cannot be shorter than one.");
	String separator;
	if (rawLength == 1) separator = "";
	else separator = "\t";
	for (int i = 0; i < size; i++) {
	    System.out.print(elementAt(i)+separator);
	    if (i%rawLength == 0) System.out.println();
	}
    }

    /**
     * Print the list elements in raws.
     * The number of elements in a raw is determined by the parameter
     **/
     public void print(int rawLength, String format) {
        if (rawLength < 1) throw new RuntimeException("Raw Length cannot be shorter than one.");
        for (int i = 0; i < size; i++) {
             System.out.printf(format,elementAt(i));
             if (i%rawLength == 0) System.out.println();
        }
	System.out.println();
     }

    public void printVerbose(int level) {
	for (int i = 0; i < size; i++)
	    if (elementAt(i) instanceof Verbose)
		System.out.println(((Verbose) elementAt(i)).verbose(level));
	    else System.out.println(elementAt(i));
    }

    public void print(MeshiWriter writer) {
	for (int i = 0; i < size; i++) {
	    writer.println(elementAt(i));
	}
	writer.flush();
    }

    public void trim() {
	Object[] newArray = new Object[size];
	for (int i = 0; i < size; i++)
	    newArray[i] = internalArray[i];
	internalArray = newArray;
	capacity = size;
    }

    /**
     * Returns the lest element of the list.
     **/
    public Object last() {
	return (elementAt(size - 1));  
    }

    /**
     * Returns modCount.
     **/
    public int modCount() {return modCount;}  
    
    public int size() { return size;}
    
    class ReversIterator implements Iterator {
        private int myModCount;
	private int current;

	public ReversIterator() {
	  current = size-1;
          myModCount = modCount;
        }

        public boolean hasNext() {
             if (myModCount != modCount)
                throw new RuntimeException("List has changed - iterator is unusable");
            return current >= 0;
													            }

        public Object next() {
            if (myModCount != modCount)
                throw new RuntimeException("List has changed - iterator is unusable");
            if (current >= 0) {
		current--;
                return internalArray[current + 1];
            }
            else return null;
         }

	public void remove() {
            throw new RuntimeException("Remove not implemented");
        }
     }

    public class MeshiListIterator implements Iterator {
	protected int myModCount;
	protected int current;
	private final Filter myFilter;
	public MeshiListIterator() {
	    current = 0;
	    myModCount = modCount;
	    myFilter = null;
	}
	public MeshiListIterator(Filter filter) {
	    current = 0;
	    myModCount = modCount;
	    this.myFilter = filter;
	}
	public boolean hasNext() {
	    if (myModCount != modCount) 
		throw new RuntimeException("\nList has changed - ModCount = "+modCount+" "+
					   "myModCount = "+myModCount+"\niterator is unusable");
	    if (myFilter != null) 
		throw new RuntimeException("\n"+"hasNext not aplicable when the iterator uses a filter");
	    return current < size;
	}
	public Object next() {
	    if (myModCount != modCount) 
		throw new RuntimeException("\nList has changed - ModCount = "+modCount+" "+
					   "myModCount = "+myModCount+"\niterator is unusable");
	    if (current < size) { 
		//		Object out = 
		current++;
		return internalArray[current - 1];
	    }
	    else return null;
	}
	public void remove() {
	    throw new RuntimeException("Remove not implemented");
	}
	public int size() { return size;}             
	public Filter listFilter() {
	    if (filter == null) throw new RuntimeException("filter is null");
	    return filter;
	}
    }  
   
    public MeshiList filter(Filter filter, MeshiList newList) {
	if (newList.size() != 0) throw new RuntimeException ("Error in MeshiList.filter(Filter filter , MeshiList newList)\n"+
							     "the new list must be empty");
	Iterator elements = iterator(); 
	Object element;
	
	newList.setComment(comment); 
	while ((element = elements.next()) != null) 
	    if (filter.accept(element)) newList.add(element);
	return newList;
    }

    /**
     * clean the the all list.
     */
    public void reset(){
    	size = 0;
    	modCount = 0;  
    }
    public boolean isEmpty() {
    	return size == 0;
    }
    static class KolDichfin implements Filter {
	public boolean accept(Object obj) {return true;}
    }

    private static class IsMeshiListFilter implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof MeshiList);
	}
    } 
}


