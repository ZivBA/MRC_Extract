package meshi.sequences;
/**
 * A container for a protein element and its number.
 **/
 public abstract class AlignmentCell {
    public final int number;
    public final Object obj;
    public final String comment;

    public AlignmentCell(Object obj, int number) {
	this.number = number;
	this.obj = obj;
	this.comment = "UNKNOWN";
    }
    
    public AlignmentCell(Object obj, int number, String comment) {
	this.number = number;
	this.obj = obj;
	this.comment = comment;
    }

    public Object object() {return obj;}
    public Object number() {return number;}
    public abstract boolean gap(); 
}
