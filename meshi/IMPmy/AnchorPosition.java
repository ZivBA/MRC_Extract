package meshi.IMPmy;

public class AnchorPosition {

	private String domainName;
	private String proteinName;
	private int resNum;
	private double x,y,z;
	private double vx,vy,vz;
	private double fx,fy,fz;
	private double snapShotX,snapShotY,snapShotZ;
	private double R;
	
	public AnchorPosition(String proteinName, String domainName, int resNum, double x, double y, double z) {
		this.domainName = domainName;
		this.proteinName = proteinName;
		this.resNum = resNum;
		this.x = x;
		this.y = y;
		this.z = z;	
		vx = 0.0;
		vy = 0.0;
		vz = 0.0;
		resetForces();
	}

	public int resNum() {
		return resNum;
	}

	public String domainName() {
		return domainName;
	}
	
	public String proteinName() {
		return proteinName;
	}
	
	public double x() {
		return x;
	}

	public void setX(double x) {
		this.x = x;
	}

	public double y() {
		return y;
	}

	public void setY(double y) {
		this.y = y;
	}

	public double z() {
		return z;
	}

	public void setZ(double z) {
		this.z = z;
	}

	public void setXYZ(double x,double y,double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	public void addX(double deltaX) {
		x+=deltaX;
	}

	public void addY(double deltaY) {
		y+=deltaY;
	}

	public void addZ(double deltaZ) {
		z+=deltaZ;
	}

	public double vx() {
		return vx;
	}

	public void setVX(double vx) {
		this.vx = vx;
	}

	public double vy() {
		return vy;
	}

	public void setVY(double vy) {
		this.vy = vy;
	}

	public double vz() {
		return vz;
	}

	public void setVZ(double vz) {
		this.vz = vz;
	}

	public void setVxyz(double vx,double vy,double vz) {
		this.vx = vx;
		this.vy = vy;
		this.vz = vz;
	}
	
	public void addVX(double deltaVX) {
		vx+=deltaVX;
	}

	public void addVY(double deltaVY) {
		vy+=deltaVY;
	}

	public void addVZ(double deltaVZ) {
		vz+=deltaVZ;
	}
	
	public void resetForces() {
		fx = 0.0;
		fy = 0.0;
		fz = 0.0;
	}
	
	public void addFX(double deltaFX) {
		fx+=deltaFX;
	}
	
	public void addFY(double deltaFY) {
		fy+=deltaFY;
	}

	public void addFZ(double deltaFZ) {
		fz+=deltaFZ;
	}
	
	public double Fx() {
		return fx;
	}
	
	public double Fy() {
		return fy;
	}
	
	public double Fz() {
		return fz;
	}

	public double R() {
		return R;
	}

	public void setR(double newR) {
		R = newR;
	}

	public void makeSnapshot() {
		snapShotX = x;
		snapShotY = y;
		snapShotZ = z;
	}
	
	public void restoreSnapshot() {
		x = snapShotX;
		y = snapShotY;
		z = snapShotZ;
	}
	
	public String toString() {
		return proteinName + "_" + domainName + "_" + resNum;
		
	}
}
