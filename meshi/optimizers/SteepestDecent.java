package meshi.optimizers;
import meshi.energy.TotalEnergy;

/**
 * This class implements a simple steepest descent minimizer, using a simple back tracking line search. 
 * The search direction is given by the gradient of the energy function, at the minimizer position.
 * The step length first tried is the previous iteration step length, multiplied by some expansion factor. 
 * If it leads to energy reduction it is choosen, otherwise it is shortened by some factor and tried again. This repeats
 * until energy reduction is achieved.
 *
 *To run this minimizer: 
 *a) Instantiate this class with the desired minimization parameters. 
 *b) Put the initial coordinates in the 'coordinates' variable at the 'energy' class.
 *c) Activate SteepestDecent.run().
 *d) Check for thrown errors to see if the minimization succeeded.
 *e) The minimized position is in the 'coordinates' variable at the 'energy' class.
 *
 * Parameters in the full constructor:
 * ----------------------------------
 * - energy - pointer to a TotalEnergy object, where the energy function is.
 * - tolerance - Minimization stops when the magnitude of the maximal gradient component drops below tolerance.
 * - maxIteration - The maximal number of iteration steps allowed
 * - reoprtEvery - The frequency of the minimization reports.
 * - initialStepLength - parameter of the line search. The first step length to be tried after the calculation of the 
 *                       first gradient. This parameter should normally be 1 unless very large gradients (such as clashhing
 *                       of VDW atoms) are expected in the first steps. In that case it should be a much smaller number. 
 * - stepSizeReduction - parameter of the line search. The step length is multiplied by this factor if no reduction 
 *                       in energy is achieved.	  
 * - stepSizeExpansion - parameter of the line search. The first step length tried is the step length from previous
 *                       line search multiplied by this factor. (Note that non-positive values to this paramater cause
 *                       special options to be called (see the SimpleStepLength class help).
 **/

public class SteepestDecent extends Minimizer {
    private SimpleStepLength lineSearch;
    private double[][] coordinates;
    private double[][] bufferCoordinates;
    private int iterationNum;
    private double lastStepLength = 1;
    private double magnitudeForce = 1000000000;
    private static final double DEFAULT_TOLERANCE = 0.00001;
    private static final int DEFAULT_MAX_ITERATION = 100000;
    private static final int DEFAULT_REPORT_EVERY = 1;
    private int reportEvery; // Default value for the frequency of the reports.
    // Default values for the line search:
    private static final double DEFAULT_INITIAL_STEP_LENGTH = 0.00000001;
    private static final double DEFAULT_STEP_SIZE_REDUCTION = 0.5;
    private static final double DEFAULT_STEP_SIZE_EXPENTION = 1.1;
    private double initialStepLength;
    private double stepSizeReduction;
    private double stepSizeExpansion;
    //Full constructor
    public SteepestDecent(TotalEnergy energy,
			  double tolerance, 
			  int maxIteration,
			  int reportEvery,
			  double initialStepLength, 
			  double stepSizeReduction, 
			  double stepSizeExpansion) {
    	super(energy,tolerance,maxIteration);
    	this.reportEvery = reportEvery;
    	this.initialStepLength = initialStepLength;
    	this.stepSizeReduction = stepSizeReduction;
    	this.stepSizeExpansion = stepSizeExpansion; 
    }
    
    // Values for the line search are taken as default in this constructor    
    public SteepestDecent(TotalEnergy energy,  
			  double tolerance, 
			  int maxIteration,
			  int reportEvery) {
	this(energy, tolerance, maxIteration, reportEvery, 
	     DEFAULT_INITIAL_STEP_LENGTH, DEFAULT_STEP_SIZE_REDUCTION, 
	     DEFAULT_STEP_SIZE_EXPENTION);
    }        

    // Default values constructor
    public SteepestDecent(TotalEnergy energy) {
	this(energy, DEFAULT_TOLERANCE, DEFAULT_MAX_ITERATION, 
	     DEFAULT_REPORT_EVERY, DEFAULT_INITIAL_STEP_LENGTH, 
	     DEFAULT_STEP_SIZE_REDUCTION, DEFAULT_STEP_SIZE_EXPENTION);
    }                
              
    private void init() {
 	coordinates = energy.coordinates();
 	lineSearch = new SimpleStepLength(energy,initialStepLength,stepSizeReduction,stepSizeExpansion);
 	bufferCoordinates = new double[coordinates.length][2];                           
 	energy.evaluate();
 	iterationNum = 0;
    }
        
    public String run() throws MinimizerException,LineSearchException {
   	init();
	isConverged = 0;
   	while ((iterationNum < maxIteration) && 
	       ((magnitudeForce =  TotalEnergy.getGradMagnitude(coordinates)) > tolerance)) {
	    for (int i = 0; i < coordinates.length; i++) {
		bufferCoordinates[i][0] = coordinates[i][0];
		bufferCoordinates[i][1] = coordinates[i][1];
	    }
	    try
		{
		    lastStepLength = lineSearch.findStepLength(bufferCoordinates);
		}
	    catch (LineSearchException e)
		{
		    if (e.code == 1)
            		throw new MinimizerException(1,"\n\n Problem in SteepestDecent."+
						     " Direction is not a descent direction. \n" + 
						     "This problem is caused by incorrect differentiation of the energy function.\n"+
						     "gradient rms is "+magnitudeForce+"\n"+
						     e.toString());
		    else
			throw e;
		}
	    iterationNum++;
	    if (iterationNum%reportEvery == 0) System.out.println(energy.report(iterationNum));
        }
        //System.out.println(energy.report(iterationNum)); //Final report  
	String finalMessage;
        if (magnitudeForce > tolerance) {
	    finalMessage = "MNC";
	    isConverged = 0;
	}
	else {
	    finalMessage = "MC";
	    isConverged = 1;
	}
	return energy.report(iterationNum)+finalMessage;
    }

    public String toString() {
	return ("SteepestDecent\n"+
		"\t maxIteration \t"+maxIteration+"\n"+
		"\t tolerance \t"+tolerance);
    }
    
    public double lastStepLength() {
    	return lastStepLength;
    }
}

		
