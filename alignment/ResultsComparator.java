package alignment;

import java.util.Comparator;

class ResultsComparator implements Comparator<Double[]>
{
	
	@Override
	public int compare(Double[] o1, Double[] o2) {
		return o2[0].compareTo(o1[0]);
	}
}

