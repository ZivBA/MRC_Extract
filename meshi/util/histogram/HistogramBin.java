package meshi.util.histogram;	
public class HistogramBin implements Comparable{
    private double centerValue;
    private int position;
    private int numberOfOccurrences;
    private double fractionOfOccurrences;
    private Histogram histogram;

    public HistogramBin(double center, int pos) {
        this(center, pos, null);
    }


    public HistogramBin(double center, int pos, Histogram histogram) {
	centerValue = center;
	position = pos;
	numberOfOccurrences = 0;
	fractionOfOccurrences = -1;
	this.histogram = histogram;
    }    
    public double centerValue() {return(centerValue);}
    public double fractionOfOccurrences(int sum) {
	fractionOfOccurrences = ((double) numberOfOccurrences)/sum;
	return(fractionOfOccurrences);
    }

    public double fractionOfOccurrences() {
        fractionOfOccurrences = ((double) numberOfOccurrences)/histogram.sum();
	return(fractionOfOccurrences);
    }

    public int numberOfOccurrences() {return(numberOfOccurrences);}	
    public HistogramBin increment() { 
	numberOfOccurrences++;
	return this;
    }
    
    public int compareTo(Object obj) {
	HistogramBin other = (HistogramBin) obj;
	if (other.fractionOfOccurrences() < fractionOfOccurrences()) return 1;
	if (other.fractionOfOccurrences() > fractionOfOccurrences()) return -1;
	return 0;
    }

    public int position() {return position;}
	
}
