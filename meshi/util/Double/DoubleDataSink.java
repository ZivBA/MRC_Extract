package meshi.util.Double;
import meshi.util.MeshiException;
import meshi.util.histogram.Histogram;

/**
 *A sink for all the junk a program may generate or collect.
 **/
public class DoubleDataSink {
    /**
     * minimum
     **/
    private double min;  

    /**
     * maximum
     **/
    private double max; 

    /**
     * sum
     **/
    private double sum;  

   /**
     * sum of squares
     **/
    private double sum2; 

    /**
     * number of data elements
     **/
    private int nElements;    

    /**
     * a histogram of the data elements
     **/
    private Histogram histogram; 

    public DoubleDataSink(Histogram histogram) {
	this();
	this.histogram = histogram;
    }

    public DoubleDataSink() {
	reset();
    }
    
    public void reset() {
	min = 10000000;
	max = -10000000;
	histogram = null;
	sum = 0;
	sum2 = 0;
	nElements = 0;
    }

	
    public void getData(double data) {
	nElements++;
	if (data < min) min = data;
	if (data > max) max = data;
	sum += data;
	sum2 += data*data;
	if (histogram != null)
	    histogram.event(data);
    }
    public int nElements() {return nElements;}
    public boolean notEmpty() {return nElements != 0;}
    public double min() { 
	if (nElements == 0) {
	    throw new MeshiException("No data elements entered the data sink");
	}
	return min;
    }
    public double max() { 
	if (nElements == 0) {
	    throw new MeshiException("No data elements entered the data sink");
	}
	return max;
    }
    public double sum() { 
	if (nElements == 0) {
	    throw new MeshiException("No data elements entered the data sink");
	}
	return sum;
    }
    public double rms() { 
	if (nElements == 0) {
	    throw new MeshiException("No data elements entered the data sink");
	}
	return Math.sqrt(sum2/nElements);
    }
    public double avg() { 
	if (nElements == 0) {
	    throw new MeshiException("No data elements entered the data sink");
	}
	return sum/nElements;
    }
    public double std() { 
	if (nElements == 0) {
	    throw new MeshiException("No data elements entered the data sink");
	}
	double avg = sum/nElements;
	double avg2 = avg*avg;
	return Math.sqrt(sum2/nElements - avg2);
    }
    public Histogram histogram() {
	return histogram;
    }
}


