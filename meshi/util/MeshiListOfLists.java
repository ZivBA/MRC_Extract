package meshi.util;
public class MeshiListOfLists {}/*extends MeshiList {
    private int elementSize;
    public MeshiListOfLists(MeshiList l1, 
			       MeshiList l2) {
	super((new Vector()).getClass());
	add(l1);
	add(l2);
    }    
    public MeshiListOfLists(MeshiList l1, 
			       MeshiList l2, 
			       MeshiList l3) {
	this(l1,l2);
	add(l3);
    }    
    public MeshiListOfLists(MeshiList l1, 
			       MeshiList l2, 
			       MeshiList l3, 
			       MeshiList l4) {
	this(l1,l2,l3);
	add(l4);
    }
    public MeshiListOfLists(MeshiList l1, 
			       MeshiList l2, 
			       MeshiList l3, 
			       MeshiList l4, 
			       MeshiList l5) {
	this(l1,l2,l3,l4);
	add(l5);
    }
    public boolean add(MeshiList list) {
	if (size() == 0)
	    {
		elementSize = list.size();
		return super.add(list);
	    }
	if (list.size() != elementSize)
	    throw new MeshiException("MeshiListOfLists error:\n Sizes must equal\n");
	return super.add(list);
    }
	
    public int elementSize() {return elementSize;}
    public MeshiList raw(int i) {
	MeshiList out = new MeshiList(); 
	MeshiList temp;
	MeshiIterator iter = Iterator();
	while ((temp = (MeshiList) iter.next()) != null)
	    out.add(temp.elementAt(i));
	return out;
    }
    public ParallelIterator parallelIterator() {
	return new ParallelIterator(this);
    }
}
		    */
