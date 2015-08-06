 package meshi.energy.inflate;
import java.util.Iterator;

import meshi.energy.EnergyElement;
import meshi.energy.Parameters;
import meshi.energy.SimpleEnergyTerm;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.AtomPair;
import meshi.optimizers.Minimizer;
import meshi.util.MeshiList;
import meshi.util.filters.Filter;

public class Inflate extends SimpleEnergyTerm{
    private static MeshiList dummy = null;
    private static Parameters parameters = new InflateParameters();
    private static AtomList atomList, atomListCopy;
    private DistanceMatrix distanceMatrix;
    private Filter filter;
    private double rmsTarget; 

    public Inflate() {}

    public Inflate(DistanceMatrix distanceMatrix, 
		   double  weight, double rmsTarget) {
	super(toArray(distanceMatrix),new InflateParametersList(), weight);
	comment = "Inflate ;)";
	this.distanceMatrix = distanceMatrix;
	on();
	this.rmsTarget = rmsTarget;
    }
    public double evaluate() {
	if (! on) return 0;
	boolean targetReached;
	try {
	    targetReached = atomList.getRms(atomListCopy) > rmsTarget;
	}
	catch (Exception e) {targetReached = true;}
	if (targetReached) Minimizer.terminator.kill("Inflate reached RMS of "+rmsTarget);
	return super.evaluate();
    }

    public void on() {
	atomList = distanceMatrix.atomList();
	atomListCopy = atomList.duplicate();
       
	createElementsList(distanceMatrix.nonBondedList());
	super.on();
    }

       public void createElementsList(DistanceList baseList) {
	elementsList = new MeshiList();

	Object baseElement;
	Parameters parameters;
	Iterator baseElements = baseList.iterator();
	while ((baseElement = baseElements.next()) != null) {
	    parameters = parametersList.parameters(baseElement);
	    if (parameters == null) handleMissingParameters(baseElement);
	    else {
		EnergyElement newElement = createElement(baseElement, parameters);
		if (! newElement.frozen())
		    elementsList.add(newElement);
	    }
	}
    }
    
    public EnergyElement createElement(Object baseElement, Parameters parameters) {
	return new InflateEnergyElement(((Distance)baseElement), distanceMatrix, weight);
    }

    public static class TargetFilter implements Filter {
	private double target;
	public TargetFilter(double target) {
	    this.target = target;
	}
	public boolean accept(Object obj) {
	    AtomPair ap = (AtomPair) obj;
	    if (ap.atom1().distanceFrom(ap.atom2()) < target) return true;
	    return false;
	}
    }
}

	
    
