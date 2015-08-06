package meshi.optimizers;
import meshi.energy.TotalEnergy;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.MeshiProgram;
import meshi.util.Terminator;

/**
 * Minimize energy according to a given set of coordinates and an energy function
 *
 **/

public abstract class Minimizer extends MeshiProgram implements KeyWords{
    public static final double DEFAULT_TOLERANCE = 0.001;
    public static final int DEFAULT_MAX_ITERATIONS = 1000;
    public static final int DEFAULT_REPORT_EVERY = 10;
    public static final int STEEPEST = 1;
    public static final int LBFGS = 2;

    protected TotalEnergy energy;
    protected double tolerance;
    protected boolean debug = MeshiProgram.debug();
    protected int reportEvery;
    protected int maxIteration;
    //Minmization terminator
    public static final Terminator terminator = new Terminator();
    /**
     * Convergence status.
     * -1 = never ran
     * 0 = last run did not converge
     * 1 = last run did converge
     **/
    protected int isConverged = -1;
 
    public Minimizer(TotalEnergy energy, double tolerance, int maxIteration) {
	this.energy = energy;
	this.tolerance = tolerance;
	this.maxIteration = maxIteration;
	terminator.reset();
    }
    public Minimizer(TotalEnergy energy, CommandList commands) {
	this.energy = energy;
	CommandList minimizerCommands = commands.firstWordFilter(MINIMIZE);
	tolerance = minimizerCommands.secondWord(TOLERANCE).thirdWordDouble();
	System.out.println("tolerance = "+tolerance);
	maxIteration = minimizerCommands.secondWord(MAX_STEPS).thirdWordInt();
	System.out.println("maxIteration = "+maxIteration);
	reportEvery = minimizerCommands.secondWord(REPORT_EVERY).thirdWordInt();
	System.out.println("reportEvery = "+reportEvery);
 	terminator.reset();
   }
    protected abstract String run() throws MinimizerException,LineSearchException;

    public String minimize() throws MinimizerException,LineSearchException {
        try{
            return run();
        }
        catch (RuntimeException re) { throw re; }
        catch (MinimizerException ex) {
            System.out.println("******************************** MINIMIZER EXCEPTION ************************\n"+
			       ex+
			       "caught in Mimimizer.minimize()");
	    ex.printStackTrace();
            energy.test();
            throw ex;
        }
        catch (LineSearchException ex) {
            System.out.println("******************************** MINIMIZER EXCEPTION ************************\n"+
			       ex+
			       "caught in Mimimizer.minimize()");
	    ex.printStackTrace();
            energy.test();
            throw ex;
        }
    }
	
    public boolean isConverged() {
	if (isConverged == -1) throw new RuntimeException(this+" never ran.");
	if (isConverged == 1) return true;
	else return false;
    }

    public TotalEnergy energy() {return energy;}
}
