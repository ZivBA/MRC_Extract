package meshi.util.crossLinking;

public class MySequence {
	
	private String title = null;
	private String seq = null;
	
	public MySequence(String title , String seq) {
		this.title = title;
		this.seq = seq;
	}
	
	public String title() {
		return title;
	}

	public String seq() {
		return seq;
	}

}
