package meshi.IMPmy;

import java.util.Vector;

public class TotalEnergy {

	private Vector<AnchorPosition> positions; // Non redundant in Anchor positions
	private DistanceConstraintList disList;
	
	public TotalEnergy(DistanceConstraintList disList) {
		this.disList = disList;
		positions = new Vector<AnchorPosition>();
		for (DistanceConstraint dc : disList) {
			if (!findPosition(dc.pos1())) {
				positions.add(dc.pos1());
			}
			if (!findPosition(dc.pos2())) {
				positions.add(dc.pos2());
			}
		}
	}
	
	public boolean findPosition(AnchorPosition pos) {
		for (AnchorPosition anotherPos : positions) {
			if (pos==anotherPos) {
				return true;
			}
		}
		return false;
	}
	
	public double evaluate(boolean withDerivatives) {
		double energy = 0.0;
		if (withDerivatives) {
			resetForces();
		}
		for (DistanceConstraint disConst : disList) {
			energy += disConst.evaluate(withDerivatives);
		}
		return energy;
	}

	public void resetForces() {
		for (AnchorPosition pos : positions) {
			pos.resetForces();
		}
	}
	
	public void snapshot() {
		for (AnchorPosition pos : positions) {
			pos.makeSnapshot();
		}
	}

	public void restoreSnapshot() {
		for (AnchorPosition pos : positions) {
			pos.restoreSnapshot();
		}
	}
	
	public void moveGrad(double gradMult) {
		for (AnchorPosition pos : positions) {
			pos.addX(gradMult*pos.Fx());
			pos.addY(gradMult*pos.Fy());
			pos.addZ(gradMult*pos.Fz());
		}
	}

	public void moveByVelocity(double dt) {
		for (AnchorPosition pos : positions) {
			pos.addX(pos.vx()*dt);
			pos.addY(pos.vy()*dt);
			pos.addZ(pos.vz()*dt);
		}
	}
	
	public void accelerate(double dt) {
		for (AnchorPosition pos : positions) {
			pos.addVX(pos.Fx()*dt);
			pos.addVY(pos.Fy()*dt);
			pos.addVZ(pos.Fz()*dt);
		}
	}

	public void scaleVelocity(double factor) {
		for (AnchorPosition pos : positions) {
			pos.setVX(pos.vx()*factor);
			pos.setVY(pos.vy()*factor);
			pos.setVZ(pos.vz()*factor);
		}
	}
	
	public Vector<AnchorPosition> positions() {
		return positions;
	}

	/**
	 * Assumes forces were calculated
	 * 
	 * returns magnitude of maximal gradient component (Inf-norm)
	 */
	public double getGradMagnitude() { 
		double maxGrad = 0.0;
		for (AnchorPosition pos : positions) {
			if (Math.abs(pos.Fx()) > maxGrad) {
				maxGrad = Math.abs(pos.Fx());
			}
			if (Math.abs(pos.Fy()) > maxGrad) {
				maxGrad = Math.abs(pos.Fy());
			}
			if (Math.abs(pos.Fz()) > maxGrad) {
				maxGrad = Math.abs(pos.Fz());
			}
		}
		return maxGrad;		
	}
	
	// Returns the root-mean-square velocity
	public double rmsVelocity() {
		double totalVelocity = 0.0;
		for (AnchorPosition pos : positions) {
			totalVelocity += (pos.vx()*pos.vx() + 
					pos.vy()*pos.vy() + 
					pos.vz()*pos.vz());
		}
		return Math.sqrt(totalVelocity/positions.size());
	}

	public void putRandomVelocities(double T) {
		for (AnchorPosition pos : positions) {
			pos.setVX(Math.random()-0.5);
			pos.setVY(Math.random()-0.5);
			pos.setVZ(Math.random()-0.5);
		}
		scaleVelocity(T/rmsVelocity());
	}
	
	public String toString() {
		String out = "";
		for (AnchorPosition pos : positions) {
			out += (pos + "   " + pos.Fx() + " " + pos.Fy() + " " + pos.Fz() + "\n");  
		}
		return out;
	}

}
