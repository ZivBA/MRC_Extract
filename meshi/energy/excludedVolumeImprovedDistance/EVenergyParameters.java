package meshi.energy.excludedVolumeImprovedDistance;
import java.util.StringTokenizer;

import meshi.energy.Parameters;
import meshi.molecularElements.Atom;
import meshi.util.filters.Filter;

public class EVenergyParameters extends Parameters {
	public final double  sigma;
	public final double  mini;
	public final double  hbEnd;
	public final boolean hb;
	public final double A,B,C,D;
	public final int smallType;
	public final int largeType;    

	public EVenergyParameters() {
		smallType = -1;
		largeType = -1;	
		sigma = -1;
		mini = -1;
		hb = false;
		hbEnd = -1;
		A = B = C = D = -1;
	}

	public EVenergyParameters(StringTokenizer st) {this(st,1.0);}

	public EVenergyParameters(StringTokenizer st,double frac) {
		int first = Atom.type(st.nextToken());
		int second = Atom.type(st.nextToken());
		if (first>second) {
			smallType = second;
			largeType = first;
		}
		else {
			smallType = first;
			largeType = second;
		}	
		int tmpInt = toInt(st.nextToken());
		if (tmpInt==0)
			hb = false;
		else
			hb = true;
		sigma = toDouble(st.nextToken())*frac;
		mini = toDouble(st.nextToken())*frac;
		hbEnd = toDouble(st.nextToken())*frac;
		
		if (mini>1e-10) { // This atom-pair has a real repulsion
			double tmp_a=0,tmp_b=0;
			double initVal = 6;
			for ( ; tmp_a>-0.01 ; initVal++) {
				tmp_a = (1/((mini-sigma)*(mini-sigma)) - initVal/((mini-sigma+0.5)*(mini-sigma+0.5)))/0.5;
				tmp_b = 1/((mini-sigma)*(mini-sigma)) - tmp_a*sigma;
			}
			A = tmp_a;
			B = tmp_b - 2*tmp_a*mini;
			C = tmp_a*mini*mini - 2*mini*tmp_b;
			D = tmp_b*mini*mini;
		}
		else
			A = B = C = D = 0.0;
	}

	public Filter isA() {return (new isA());}

	public Parameters create(StringTokenizer stringTokenizer) {
		return (new  EVenergyParameters(stringTokenizer));
	}

	private  class isA implements Filter {
		public boolean accept(Object obj) {
			return (obj instanceof EVenergyParameters);
		}
	}

	public String toString() {return ""+sigma;}
}
