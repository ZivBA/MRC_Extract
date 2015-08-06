package meshi.util.mathTools;

/**
 * This class calculates sigma for a set of parameters. It gives a derivable sigmoid function
 * of the following form:
 *
 * x<0   - not defined
 * 0<x<start - 1.0
 * start<x<p1 - a quadratic decent from 1.0 to the value of valAtp1
 * p1<x<p2 - a cubic spline decent from (p1,valAtp1) to (p2,valAtp2)
 * p2<x<end - a quadratic decent from (p2,valAtp2) to (end,0.0).  
 * end<x - 0.0
 * 
 * Use s(x), and s_tag(x) to get 's' and 's_tag' values for a certain x. 
 *
 * Note: The sigmoid function is derivable for all x=<0.
 **/   

public class HB_Sigma {	
	// Sigma parameters
	private double start,p1,p2,end,valAtp1,valAtp2;
	
	// Auxillary variables:
	// --------------------
	// quadratic at [p2,end]
	private double partC_a;
	
	// Spline section
	private double deriv_at_p2,deriv_at_p1;
	private double a,b,c,d,xMinusP1;

	// quadratic at [start,p1]
	private double partA_a;

	public HB_Sigma(double start , double p1, double p2, double end,
								double valAtp1, double valAtp2) {
		this.start = start;
		this.p1 = p1;
		this.p2 = p2;
		this.end = end;
		this.valAtp1 = valAtp1;
		this.valAtp2 = valAtp2;
		
		// quadratic at [p2,end]
		partC_a = valAtp2/((p2-end)*(p2-end));
		deriv_at_p2 = 2.0*partC_a*(p2-end);
		
		// quadratic at [start,p1]
		partA_a = (1.0-valAtp1)/((p1-start)*(p1-start));
		deriv_at_p1 = -2.0*partA_a*(p1-start);
		
		// spline
		d = valAtp1;
		c = deriv_at_p1;
		// for a,b I solve a 2 variable linear equation
		double A = (p2-p1)*(p2-p1)*(p2-p1);
		double B = (p2-p1)*(p2-p1);
		double C = valAtp2-d-c*(p2-p1);
		double D = 3*(p2-p1)*(p2-p1);
		double E = 2*(p2-p1);
		double F = deriv_at_p2 - c;
		a = (C*E-B*F)/(A*E-B*D);
		b = (A*F-D*C)/(A*E-B*D);
	}
	
	/**
	 *  The value of sigma
	 **/
	public double s(double x) {
    	if (x>=end) {
    		return 0.0;
    	}
    	if (x>=p2) {
    		return partC_a*(x-end)*(x-end);
    	}
    	if (x>=p1) {
    		xMinusP1 = x-p1;
    		return a*xMinusP1*xMinusP1*xMinusP1 + b*xMinusP1*xMinusP1 + c*xMinusP1 + d;
    	}
    	if (x>start) {
    		return 1.0-partA_a*(x-start)*(x-start);
    	}
    	if (x>=-1e-12) {
    		return 1.0;
    	}
    	throw new RuntimeException("A negative x-value is not possible: " + x);  
	}

	/**
	 *  The derivative of sigma
	 **/
	public double s_tag(double x) {
    	if (x>=end) {
    		return 0.0;
    	}
    	if (x>=p2) {
    		return 2.0*partC_a*(x-end);
    	}
    	if (x>=p1) {
    		xMinusP1 = x-p1;
    		return 3.0*a*xMinusP1*xMinusP1 + 2.0*b*xMinusP1 + c;
    	}
    	if (x>start) {
    		return (-2.0)*partA_a*(x-start);
    	}
    	if (x>=-1e-12) {
    		return 0.0;
    	}
    	throw new RuntimeException("A negative x-value is not possible: " + x);  
	}

}