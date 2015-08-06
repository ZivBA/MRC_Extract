package meshi.util.histogram;

import java.util.Arrays;
import java.util.Iterator;

import meshi.util.MeshiIterator;
import meshi.util.SortableMeshiList;
import meshi.util.filters.Filter;


public class Histogram extends SortableMeshiList{
    private double lowBound;
    private double highBound;
    private int nBins;
    private double histFactor;
    private int belowLow,aboveHigh;
    private int sum;
    public Histogram(double low, double high,int n) {
	super(new IsBin());
	lowBound = low;
	highBound =  high;
	nBins = n;
	if (lowBound >= highBound ) 
	    throw new RuntimeException("Histogram constructor Error:\n"+
				       "lowBound "+lowBound+" >= highBound "+highBound+" \n");
	if (nBins <= 0) new RuntimeException("Histogram constructor Error:\n"+
				       "nBins = "+nBins+" \n");
	belowLow = aboveHigh = 0;
	histFactor = nBins/(highBound - lowBound);
	sum = 0;
	for (int i = 0; i < nBins; i++)
	    add(new HistogramBin(centerOfBin(i),i,this));
    }
    private int pos(double x) {
	if (lastSorted != -1) 
	    throw new RuntimeException("Histogram.pos inapplicable after sort"); 
	if (x >= highBound) return -2;
	if (x < lowBound) return -1;
	return((new Double((x-lowBound)*histFactor)).intValue());
    }
    private double centerOfBin(int bin) {
	if (lastSorted != -1) 
	    throw new RuntimeException("Histogram.centerOfBin inapplicable after sort"); 
        return((bin+0.5)/histFactor+lowBound);
    }
			  
    public void event(double x) {
	if (lastSorted != -1) 
	    throw new RuntimeException("Histogram.event inapplicable after sort"); 
	int bin = pos(x);
	if (bin == -2) aboveHigh++;
	else 
	    if (bin  == -1) belowLow++;
	    else {
		if (bin == size()) // may result from numerical instability (10E-16)+10 = 10
		    bin = size()-1;
		set(bin,((HistogramBin) elementAt(bin)).increment()); 
		/* set is used in the else term above in order to use the  
		   concurrentModificationException mechanizm. */
		sum++;
		}
    }
    public HistogramBin HistogramBinAt(int pos){
	if (lastSorted != -1) 
	    throw new RuntimeException("Histogram.HistogramBinAt inapplicable after sort"); 
	return (HistogramBin) elementAt(pos);
    }
				      
    public HistogramBin HistogramBinAt(double x) {
	if (lastSorted != -1) 
	    throw new RuntimeException("Histogram.HistogramBinAt inapplicable after sort"); 
	return(HistogramBinAt(pos(x)));
    }

    public double lowBound() {return(lowBound);}
    public double highBound() {return(highBound);}
    public int nBins() {return(nBins);}
    public HistogramBin highestPeak() {
	int max = -100;
	int pos = -1;
	int value;
	Iterator iter = iterator();
	HistogramBin temp,maxBin = new HistogramBin(-1,-1);
	while ((temp = (HistogramBin) iter.next()) != null)
	    {
		value = temp.numberOfOccurrences();
		if (max < value)
		    {
			max = value;
			maxBin = temp;
		    }
	    }
	return maxBin;
    }
    
    static class BinCreator {
	public Object create() {return new HistogramBin(0.0, 0);}
    }
    static class IsBin implements Filter {
	public boolean accept(Object obj) {
		return (obj instanceof HistogramBin);
       }
    }
    public void print() {
	HistogramBin bin;
	Iterator iter = iterator();
	while ((bin = (HistogramBin) iter.next()) != null)
	    System.out.println(bin.centerValue()+"\t"+
			       bin.fractionOfOccurrences());
    }
    public boolean sortable() {return true;}
    public double sum() {return sum;}
    
    public Iterator representativesIterator(double fractionThreshold) {
	return new RepresentativesIterator(fractionThreshold);
    }
    public HistogramBin representative(double fraction) {
	    Iterator bins = iterator();
	    HistogramBin bin,prev = null;
	    
	    double sum = 0;
	    while (((bin = (HistogramBin) bins.next()) != null) & 
		   (sum < fraction)) {
		prev = bin;
		sum +=bin.fractionOfOccurrences();
	    }
	    return prev;
    }
    private  class RepresentativesIterator implements Iterator {
	MeshiIterator internal;
	ThresholdFilter filter;
	public RepresentativesIterator(double fractionThreshold) {
	    double[] values = new double[size];
	    Iterator bins = iterator();
	    HistogramBin bin;
	    
	    int i = 0;
	    while ((bin = (HistogramBin) bins.next()) != null) {
		values[i] = bin.fractionOfOccurrences();
		i++;
	    }
	    Arrays.sort(values);
	    
	    double sum = 0;
	    for (i = size-1; (i>= 0) & sum <= fractionThreshold; i--) {
		sum += values[i];
	    }
	    filter = new ThresholdFilter(values[i]);

	    internal = meshiIterator();
	}
	
	public Object next() {
	    return internal.next(filter);
	}

	public boolean hasNext() {
	    return internal.hasNext();
	}
	public void remove() {
	    internal.remove();
	}
	
	private class ThresholdFilter implements Filter {
	    double threshold; 
	    public ThresholdFilter(double threshold) {
		this.threshold = threshold;
	    }
	    public boolean accept(Object obj) {
		HistogramBin bin = (HistogramBin) obj;
		return bin.fractionOfOccurrences()>threshold;
	    }
	}
    }
}
