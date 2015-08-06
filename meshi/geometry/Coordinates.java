package meshi.geometry;

public class Coordinates{
    protected double[] x,y,z;
    public  Coordinates() {
	this(0.0, 0.0, 0.0);
    }
    public  Coordinates(double x, double y, double z) {
	this.x = new double[2];
	this.y = new double[2];
	this.z = new double[2];
	this.x[0] = x;	
	this.y[0] = y;	
	this.z[0] = z;
	this.x[1] = 0.0;
	this.y[1] = 0.0;
	this.z[1] = 0.0;
    }
    public  Coordinates(Coordinates coordinates) { 
	this(coordinates.x(),
	     coordinates.y(),
	     coordinates.z());
    }
    public final double x() {return x[0];}
    public final double y() {return y[0];}
    public final double z() {return z[0];}

    public final double fx() {return x[1];}
    public final double fy() {return y[1];}
    public final double fz() {return z[1];}

    public void setX(double x) {this.x[0] = x;}
    public void setY(double y) {this.y[0] = y;}
    public void setZ(double z) {this.z[0] = z;}

    public void setFx(double fx) {this.x[1] = fx;}
    public void setFy(double fy) {this.y[1] = fy;}
    public void setFz(double fz) {this.z[1] = fz;}

    public void addToX(double addMe) {x[0] += addMe;}
    public void addToY(double addMe) {y[0] += addMe;}
    public void addToZ(double addMe) {z[0] += addMe;}

    public void addToFx(double addMe) {x[1] += addMe;}
    public void addToFy(double addMe) {y[1] += addMe;}
    public void addToFz(double addMe) {z[1] += addMe;}

    public void resetForces() {
	x[1] = 0.0;
	y[1] = 0.0;
	z[1] = 0.0;
    }
    public final double[] X() {return x;}
    public final double[] Y() {return y;}
    public final double[] Z() {return z;}

    public boolean equals(Object obj) {
	Coordinates other = (Coordinates) obj;
	return
            ((this.z() == other.z())&
	     (this.y() == other.y())&
	     (this.x() == other.x()));
    }
}
