package alignment;

import java.util.Vector;

public abstract class Sequence extends Vector<Position> {

	public int length() {
		return size();
	}
}
