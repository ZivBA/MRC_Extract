package meshi.optimizers;
import meshi.energy.TotalEnergy;

/**
 * This class implements a nonlinear Conjugate Gradient minimizer 
 * 
 * PR+ algorithm is implemented.
 * 
 * Fletcher-Reeves algorithm is the base (FR-CG 5.4 in book below) 
 * using betaPR (PR-CG algorithm) (Polak-Ribiere) 5.43 and using 5.44 (non-negative beta: PR+) pp. 122.
 * with posible restarts every n steps.
 * according to the scheme in: Numerical Optimization by J. Nocendal & 
 * S. J. Wright, Springer 1999, pp 120-122.
 *
 **/

public class ConjugateGradient extends Minimizer {
    protected LineSearch lineSearch;
    protected double[][] coordinates;
    protected double[][] bufferCoordinates;
    protected double[] P; // search direction
    protected double[] G; // The (-) gradients at iteration k	
    protected double beta; // beta at iteration k+1
    
	private static final int DEFAULT_RESTART_EVERY = 0;
	protected int restartEvery;

    private int iterationNum;
    private double magnitudeForce;
    
    private static final double DEFAULT_TOLERANCE = 0.00001;
    private static final int DEFAULT_MAX_ITERATION = 100000;
    private static final int DEFAULT_REPORT_EVERY = 10;
    private static  double initStepSteepestDecent = 0.0001;
    private static  double stepSizeReductionSteepestDecent = 0.5;
    private static  double stepSizeExpansionSteepestDecent = 2.1;
    private static  int numStepsSteepestDecent = 100;
    private static  int maxSteepestDecent = 6;
    private int nSteepestDecent = 0;

    private int reportEvery; // Default value for the frequency of the reports.
    
	// Wolf conditions line search parameters
	private double c1;
	private double c2;
	private double extendAlphaFactorWolfSearch;
	private int maxNumEvaluationsWolfSearch;
	private static final double DEFAULT_C1 = 1e-3;
	private static final double DEFAULT_C2 = 0.4;
	private static final double DEFAULT_EXTENDED_ALPHA_FACTOR_WOLF_SEARCH = 3.0;
	private static final int DEFAULT_MAX_NUM_EVALUATIONS_WOLF_SEARCH = 10;

    //Full constructor
	public ConjugateGradient(TotalEnergy energy,
				 double tolerance,
				 int maxIteration,
				 int reportEvery,
				 double c1,
				 double c2,
				 double extendAlphaFactorWolfSearch,
				 int maxNumEvaluationsWolfSearch,
				 int restartEvery) {
		super(energy, tolerance, maxIteration);
		this.reportEvery = reportEvery;
		this.c1 = c1;
		this.c2 = c2;
		this.extendAlphaFactorWolfSearch = extendAlphaFactorWolfSearch;
		this.maxNumEvaluationsWolfSearch = maxNumEvaluationsWolfSearch;
		this.restartEvery = restartEvery;
	}
    
    // Values for the line search are taken as default in this constructor    
    public ConjugateGradient(TotalEnergy energy,  
			  double tolerance, 
			  int maxIteration,
			  int reportEvery) {
	this(energy, tolerance, maxIteration, reportEvery, 
			DEFAULT_C1, DEFAULT_C2,
			DEFAULT_EXTENDED_ALPHA_FACTOR_WOLF_SEARCH,
			DEFAULT_MAX_NUM_EVALUATIONS_WOLF_SEARCH,
			DEFAULT_RESTART_EVERY);
    }        

	// another constructor to specify the restart
	public ConjugateGradient(TotalEnergy energy,  
			  double tolerance, 
			  int maxIteration,
			  int reportEvery,
			  int restartEvery) {
	this(energy, tolerance, maxIteration, reportEvery, 
			DEFAULT_C1, DEFAULT_C2,
			DEFAULT_EXTENDED_ALPHA_FACTOR_WOLF_SEARCH,
			DEFAULT_MAX_NUM_EVALUATIONS_WOLF_SEARCH,
			restartEvery);
	}        

    // Default values constructor
    public ConjugateGradient(TotalEnergy energy) {
	this(energy, DEFAULT_TOLERANCE, DEFAULT_MAX_ITERATION, 
	     DEFAULT_REPORT_EVERY, DEFAULT_C1, DEFAULT_C2,
		 DEFAULT_EXTENDED_ALPHA_FACTOR_WOLF_SEARCH,
		 DEFAULT_MAX_NUM_EVALUATIONS_WOLF_SEARCH,
		 DEFAULT_RESTART_EVERY);
    }                
              
	private void init() {
		coordinates = energy.coordinates();
		lineSearch =
			new WolfConditionLineSearch(energy,
										c1,
										c2,
										extendAlphaFactorWolfSearch,
										maxNumEvaluationsWolfSearch);
		bufferCoordinates = new double[coordinates.length][2];
		P = new double[coordinates.length];  
		G = new double[coordinates.length];

		energy.evaluate();
		iterationNum = 0;
		
		//init P and G
		for (int i = 0; i < coordinates.length; i++) { 
			P[i] = coordinates[i][1];				
			G[i] = coordinates[i][1];
		}			
	}
        
	public String run() throws MinimizerException, LineSearchException{
	    
		init();
		
		//main loop
		while ((iterationNum < maxIteration) && 
		       ((magnitudeForce = TotalEnergy.getGradMagnitude(coordinates)) > tolerance)) {

			// do line search
			for (int i = 0; i < coordinates.length; i++) {
				bufferCoordinates[i][0] = coordinates[i][0];
				bufferCoordinates[i][1] = P[i];
			}
			try {
			    lineSearch.findStepLength(bufferCoordinates);
			}
			catch (LineSearchException lse) {
			    if (nSteepestDecent >= maxSteepestDecent)
				throw new MinimizerException(10,"too many kickstarts "+lse);
			    SteepestDecent steepestDecent = new SteepestDecent(energy,tolerance,numStepsSteepestDecent,
									       reportEvery,initStepSteepestDecent,
									       stepSizeReductionSteepestDecent,
									       stepSizeExpansionSteepestDecent);
			    try {
				if (iterationNum != 0) 
				    System.out.println("# A kick start has occurred in iteration:"+iterationNum);
				steepestDecent.run();
				nSteepestDecent++;
			    }
			    catch(MinimizerException me)
				{
				    if (me.code != 2)
					throw me; //Differentiation error, or some other irrecoverable error 
				}
			}

			// calculate beta
			double betaSum = 0, normalSumSquares = 0;
			for (int i = 0; i < coordinates.length; i++) {
				betaSum += coordinates[i][1] * (coordinates[i][1] - G[i]);
				normalSumSquares += G[i]*G[i];
			}				
			beta = betaSum / normalSumSquares; 
			if (beta < 0 || (restartEvery!=0 && iterationNum!=0 && iterationNum%restartEvery == 0)) {
				//System.out.println("beta reset at "+iterationNum+" (was "+beta+")");
				beta = 0;
			}
			
			// calculate Pk+1, G
			for (int i = 0; i < coordinates.length; i++) { 			
				P[i] = coordinates[i][1] + beta * P[i];
				G[i] = coordinates[i][1];
			}
				
			iterationNum++;
			if (iterationNum % reportEvery == 0)
				System.out.println(energy.report(iterationNum));
		}

		//		System.out.println(energy.report(iterationNum)); //Final report     
		return energy.report(iterationNum);
	}
    
    public String toString() {return "ConjugateGradient";}
}
	
		
