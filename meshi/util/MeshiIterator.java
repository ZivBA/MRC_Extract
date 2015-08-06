package meshi.util;

import java.util.ConcurrentModificationException;
import java.util.Iterator;

import meshi.util.filters.Filter;

/**
 * An iterator for nested loops.
 * Note that for simple loops the Iterator returned by {@link meshi.util.MeshiList#iterator} is
 * a better choice as it runs faster.
 * <p>
 * Conceptualy each MeshiIterator has an internal MeshiIterator over the same list.
 * In practice of course nested iterators are instantiated on demand.
 * A nested iterator that looped over all list elements dies, and is replaced by a new iterator in the 
 * next step of the nesting iterator.
 * </p><p>
 * <b>Example:</b><br>
 * let <b>ml</b> be a list with the Integer Objects 1-4<br>
 * <pre>
 * 	  MeshiIterator mi = ml.meshiIterator();
 *	  while((temp = mi.next()) != null)
 * 	       while ((temp1 = mi.nestedNextFrom()) != null)
 *                    while ((temp2 = mi.nested().nestedNextFrom()) != null)
 * 			   System.out.println(temp+" "+temp1+" "+temp2);
 *</pre>
 * will result in the following output:<br>
 * 1 2 3<br>
 * 1 2 4<br>
 * 1 3 4<br>
 * 2 3 4<br>
 **/
public class MeshiIterator implements Iterator{
    /**
     * The list over which we iterate.  
     **/
    protected MeshiList list; 

    /**
     * of the list.
     **/
    private int size; // of the list

    /**
     * of the list at the time of the
     * iterator creation. The list modCount
     * is not expected to change while an iterator
     * is living.
     **/
    private int modCount;

    /**
     * of the iterator.
     **/
    private int currentPosition; 

    /**
     * the iterator iterates between <b>first</b> and <b>last</b> elements of the list.
     **/ 
    private int first; 

    /**
     * the iterator iterates between <b>first</b> and <b>last</b> elements of the list.
     **/ 
    private int last;
    private boolean dead; // A flag which mean that farther iterations are
                          // not needed. 
    /**
     * An iterator over the same list. Used for nested loops.
     **/
    private MeshiIterator nested; 
 
    /** 
     *Instantiates a MeshiIterator over the <b>list</b>.
     **/
    public MeshiIterator(MeshiList list) { 
	this.list = list;                    
	size = list.size();
	first = 0;
	last = size -1;
	modCount = list.modCount();
	currentPosition = -1;
	dead = false;
	nested = null;
    }

    //methods for easy houskeeping
    /**
     * that the list has not been modified since this iterator was insantiated
     **/
    private void check() { // that the list has not been modified
	if (list == null) return; 
	if (modCount != list.modCount()) 
	    MeshiIteratorError.ConcurrentModification(list);
    }

    /**
     *  Starting point for iterations.
     */
    protected  void setFirst(int first) { 
	check();
	this.first = first;
	currentPosition = first-1;
	if (first < 0) dead = true; 
    }

    /**
     *  Last point for iterations.
     */
    protected  void setLast(int last) { // end point
	check();
	this.last = last;
	if (last >= size ) dead = true; 
    }
    // end of housekeeping

    public boolean hasNext() {
	check();
	if (dead) return false;
	if (currentPosition >= last) return false;
	return true;
    }

    private int currentPosition() {
	check();
	if (dead || (currentPosition >= last)) return -1;
	return currentPosition;
    }

    /**
     * Next object of the list.
     **/
    public Object next() {
	currentPosition ++;
	return current();
    }

    /**
     * :) .
     **/
    private Object next1() {
	currentPosition ++;
	return current();
    }

    /**
     * The object in position previous+step.
     **/
    public Object next(int step) {
	currentPosition += step;
	return current();
    }

    /**
     * Next object that pass the filter.
     **/
    public Object next(Filter filter) {
	Object nextObject;
	while ((nextObject = next1()) != null)
	    if (filter.accept(nextObject)) return nextObject;
	return null; 
    }

    /**
     * The object at the current position.
     **/
    private Object current()  {
	check();
	if (dead) return null;
	if (currentPosition > last) {
	    kill();
	    return null;
	    }
	return list.elementAt(currentPosition);
    }

    /**
     * the iterator (hasNext() returns false)
     **/
    public void kill() {
	check();
	dead = true;
    }

    public boolean dead() { return dead;}
    
    /**
     * Returns a nested iterator.
     **/
     public  MeshiIterator nested() {
 	return nested;
     }
    
    /**
     * Returns the next object of the nested iterator that pass the <b>filter</b>.
     * The nested iterator iterates between position <b>from</b>of the list to
     * position <b>last</b>.<br>
     * Note that this method also instantiate the nested iterator in each iteration of the 
     * nesting iterator.
     **/
    private Object nestedNext(int from, int to,Filter filter) {
	check();
	if ((nested == null) || nested.dead()) {
	    nested = (MeshiIterator) list.meshiIterator();
	    nested.setFirst(from); 
	    nested.setLast(to);  
	}
	return nested.next(filter);
    }

    /**
     * Returns the next object of the nested iterator.
     * The nested iterator iterates between position <b>from</b>of the list to
     * position <b>last</b>.<br>
     * Note that this method also instantiate the nested iterator in each iteration of the 
     * nesting iterator.
     **/
    private Object nestedNext(int from, int to) {
	check();
	if ((nested == null) || nested.dead()) {
	    nested = (MeshiIterator) list.meshiIterator();
	    nested.setFirst(from); 
	    nested.setLast(to);  
	}
	return nested.next();
    }
    
    /**
     * Returns the next object of the nested iterator.
     * The nested iterator starts from the current position of the nesting iterator. 
     * That is it behaves like the nested loop in:<br>
     *                   for(i = 0; i < n; i++)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    <= nesting loop<br>
     * &nbsp;&nbsp;                  for(j = i + 1; j < n; j++)                  <= nested loop<br>
     * Note that this method also instantiate the nested iterator in each iteration of the 
     * nesting iterator.
     **/
    public Object nestedNextFrom() {
	return nestedNext(currentPosition+1,last);
    }

    /**
     * Returns the next object of the nested iterator.
     * The nested iterator dies when it reaches one position before the current position of 
     * the nesting iterator. 
     * That is it behaves like the nested loop in:<br>
     *                   for(i = 0; i < n; i++)&nbsp;&nbsp;    <= nesting loop<br>
     * &nbsp;&nbsp;                  for(j = 0; j < i; j++)    <= nested loop<br>
     * Note that this method also instantiate the nested iterator in each iteration of the 
     * nesting iterator.
     **/
    public Object nestedNextTo() {
	return nestedNext(first,currentPosition-1);
    }

    /**
     * Returns the next object of the nested iterator.
     * The nested iterator dies when it reaches the current position of the nesting iterator. 
     * That is it behaves like the nested loop in:<br>
     *                   for(i = 0; i < n; i++)&nbsp;&nbsp;    <= nesting loop<br>
     * &nbsp;&nbsp;                  for(j = 0; j <= i; j++)   <= nested loop<br>
     * Note that this method also instantiate the nested iterator in each iteration of the 
     * nesting iterator.
     **/
    public Object nestedNextToInclude() {
	return nestedNext(first,currentPosition);
    }

    /**
     * Returns the next object of the nested iterator that pass the filter.
     * Note that this method also instantiate the nested iterator in each iteration of the 
     * nesting iterator.
     **/
    public Object nestedNext(Filter filter) {
	return nestedNext(first,last,filter);
    }

     /**
     * Returns the next object of the nested iterator.
     **/
    public Object nestedNext() {
	return nestedNext(first,last);
    }

    /**
     * Returns the next object of the nested iterator that pass the filter.
     * The nested iterator starts from the current position of the nesting iterator. 
     * That is it behaves like the nested loop in:<br>
     *                   for(i = 0; i < n; i++)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    <= nesting loop<br>
     * &nbsp;&nbsp;                  for(j = i + 1; j < n; j++)                  <= nested loop<br>
     * Note that this method also instantiate the nested iterator in each iteration of the 
     * nesting iterator.
     **/
    public Object nestedNextFrom(Filter filter) {
	return nestedNext(currentPosition+1,last,filter);
    }

    /**
     * Returns the next object of the nested iterator that pass the filter.
     * The nested iterator dies when it reaches the current position of the nesting iterator. 
     * That is it behaves like the nested loop in:<br>
     *                   for(i = 0; i < n; i++)&nbsp;&nbsp;    <= nesting loop<br>
     * &nbsp;&nbsp;                  for(j = 0; j < i; j++)    <= nested loop<br>
     * Note that this method also instantiate the nested iterator in each iteration of the 
     * nesting iterator.
     **/
    public Object nestedNextTo(Filter filter) {
	return nestedNext(first,currentPosition-1,filter);
    }

    /**
     * Not implemented.
     **/
    public void remove() {
      throw new MeshiException("remove not implemented yet in Meshi iterator");
    }
}
class MeshiIteratorError {
    public  static void ConcurrentModification(MeshiList list) throws ConcurrentModificationException {
	throw new ConcurrentModificationException("MeshiIterator error\n"+
						   "modified list :"+list.toString()+"\n");
    }
}
    
