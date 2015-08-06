package meshi.optimizers;

public class LineSearchException extends Exception {
	
	public int code;
	
	public LineSearchException(int code, String msg) {
		super(msg);		
		this.code = code;
	}
}



