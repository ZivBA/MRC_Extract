package meshi.IMPmy;

public class DistanceConstraint {

	private AnchorPosition pos1;
	private AnchorPosition pos2;
	private double weight;
	private double targetDis = Double.NEGATIVE_INFINITY;
	public final DistanceConstraintType constraintType;

	// For calculations
	private double energy;
	private double dEdD;
	private double dx;	
	private double dy;	
	private double dz;
	private double dis;
	private double invDistance;
	private double dDistanceDx;
	private double dDistanceDy;
	private double dDistanceDz;
	private double rEV = 0.0;

	public DistanceConstraint(AnchorPosition pos1 , AnchorPosition pos2, DistanceConstraintType constraintType, double weight) {
		this.pos1 = pos1;
		this.pos2 = pos2;
		this.constraintType = constraintType;
		this.weight = weight;
		if (constraintType == DistanceConstraintType.RIGID_BODY) {
			targetDis = Math.sqrt(Math.pow((pos1.x()-pos2.x()),2.0) + Math.pow((pos1.y()-pos2.y()),2.0) + Math.pow((pos1.z()-pos2.z()),2.0));
		}
		rEV = (pos1.R() + pos2.R());
	}
		
	
	// This code for parsed sequences.
	public double evaluate(boolean withDerivatives) {
		double CONNECTIVITY_MAX = 12.1; // Angs  // For parding of 25 res I used conn_max of 13.0 and rEV of 12.0
		double XL_MAX = 20.0; // Angs
		rEV = 12.0;
		energy = 0.0;
		dx  = pos2.x()-pos1.x(); 	
		dy  = pos2.y()-pos1.y();	
		dz  = pos2.z()-pos1.z();
		dis = Math.sqrt(dx*dx + dy*dy + dz*dz);
		if (withDerivatives) {
			invDistance = 1.0/dis;
			dDistanceDx = dx*invDistance;
			dDistanceDy = dy*invDistance;
			dDistanceDz = dz*invDistance;
		}
		
		switch (constraintType) {
		case CONNECTIVITY:	
			if (dis<CONNECTIVITY_MAX) {
				energy = dEdD = 0.0;
			} 
			else { 
				double alpha = 1.0;
				double inv1 = 1.0/(1 + (dis-CONNECTIVITY_MAX)*alpha);
				energy = weight*alpha*(dis-CONNECTIVITY_MAX)*(dis-CONNECTIVITY_MAX)*inv1;
				dEdD = weight*alpha*(dis-CONNECTIVITY_MAX)*(2*inv1 - alpha*(dis-CONNECTIVITY_MAX)*inv1*inv1);
			}
			break;

		case CROSS_LINK: 
			if (dis<XL_MAX) {
				energy = dEdD = 0.0;
			} 
			else { 
				double alpha = 1.0;
				double inv1 = 1.0/(1 + (dis-XL_MAX)*alpha);
				energy = weight*alpha*(dis-XL_MAX)*(dis-XL_MAX)*inv1;
				dEdD = weight*alpha*(dis-XL_MAX)*(2*inv1 - alpha*(dis-XL_MAX)*inv1*inv1);
			}
			break;

		case RIGID_BODY: 
			energy = weight*(dis-targetDis)*(dis-targetDis);
			dEdD = 2*weight*(dis-targetDis);
			break;

		case EXCLUDED_VOLUME:
			if (dis>rEV) {
				energy = dEdD = 0.0;
			} 
			else { 
				energy = weight*(dis-rEV)*(dis-rEV);
				dEdD = 2.0*weight*(dis-rEV);
			}
			break;
			
		default:	 
			energy = dEdD = 0.0;
			break;
		}
		
		if (withDerivatives) {
			pos1.addFX(dEdD*dDistanceDx); // These are forces not derivatives
			pos1.addFY(dEdD*dDistanceDy);
			pos1.addFZ(dEdD*dDistanceDz);
			pos2.addFX(-dEdD*dDistanceDx);
			pos2.addFY(-dEdD*dDistanceDy);
			pos2.addFZ(-dEdD*dDistanceDz);
		}
		
		return energy;
	}
	

//	// This code for structured domain list.
//	public double evaluate(boolean withDerivatives) {
//		energy = 0.0;
//		dx  = pos2.x()-pos1.x(); 	
//		dy  = pos2.y()-pos1.y();	
//		dz  = pos2.z()-pos1.z();
//		dis = Math.sqrt(dx*dx + dy*dy + dz*dz);
//		if (withDerivatives) {
//			invDistance = 1.0/dis;
//			dDistanceDx = dx*invDistance;
//			dDistanceDy = dy*invDistance;
//			dDistanceDz = dz*invDistance;
//		}
//		
//		switch (constraintType) {
//		case CONNECTIVITY:	
//			if (dis<5.0) {
//				energy = dEdD = 0.0;
//			} 
//			else { 
//				double alpha = 1.0;
//				double inv1 = 1.0/(1 + (dis-5.0)*alpha);
//				energy = weight*alpha*(dis-5.0)*(dis-5.0)*inv1;
//				dEdD = weight*alpha*(dis-5.0)*(2*inv1 - alpha*(dis-5.0)*inv1*inv1);
//			}
//			break;
//
//		case CROSS_LINK: 
//			if (dis<30.0) {
//				energy = dEdD = 0.0;
//			} 
//			else { 
//				double alpha = 1.0;
//				double inv1 = 1.0/(1 + (dis-30.0)*alpha);
//				energy = weight*alpha*(dis-30.0)*(dis-30.0)*inv1;
//				dEdD = weight*alpha*(dis-30.0)*(2*inv1 - alpha*(dis-30.0)*inv1*inv1);
//			}
//			break;
//
//		case RIGID_BODY: 
//			energy = weight*(dis-targetDis)*(dis-targetDis);
//			dEdD = 2*weight*(dis-targetDis);
//			break;
//
//		case EXCLUDED_VOLUME:
//			if (dis>rEV) {
//				energy = dEdD = 0.0;
//			} 
//			else { 
//				energy = weight*(dis-rEV)*(dis-rEV);
//				dEdD = 2.0*weight*(dis-rEV);
//			}
//			break;
//			
//		default:	 
//			energy = dEdD = 0.0;
//			break;
//		}
//		
//		if (withDerivatives) {
//			pos1.addFX(dEdD*dDistanceDx); // These are forces not derivatives
//			pos1.addFY(dEdD*dDistanceDy);
//			pos1.addFZ(dEdD*dDistanceDz);
//			pos2.addFX(-dEdD*dDistanceDx);
//			pos2.addFY(-dEdD*dDistanceDy);
//			pos2.addFZ(-dEdD*dDistanceDz);
//		}
//		
//		return energy;
//	}
	
	
	public void setWeight(double newWeight) {
		weight = newWeight;
	}
	
	public AnchorPosition pos1() {
		return pos1;
	}

	public AnchorPosition pos2() {
		return pos2;
	}

	public String toString() {
		return pos1.proteinName()+"_"+pos1.domainName()+"_"+pos1.resNum()+" <----    "+constraintType+"    ----> "+pos2.proteinName()+"_"+pos2.domainName()+"_"+pos2.resNum() + "                Energy: " + evaluate(false); 
	}

}
