package meshi.geometry;

public class FrozenCoordinates extends Coordinates{
    private final double INFINITY = 1/0.0;
    public   FrozenCoordinates(Coordinates coordinates) { 
	super(coordinates.x(),
	      coordinates.y(),
	      coordinates.z());
    	x[1] = INFINITY;
	y[1] = INFINITY;
	z[1] = INFINITY;
    }
    public void setFx(double fx) {}
    public void setFy(double fy) {}
    public void setFz(double fz) {}
    public void resetForces() {}
}
