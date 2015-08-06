package meshi.optimizers;

public class MinimizerException extends Exception {
	
	public final int code;
	
	public MinimizerException(int code, String msg) {
		super(msg);		
		this.code = code;
	}
}

// *********************************************************************************************************
/* The following code could be useful to test if the exception was thrown because of diffrentiation problem:
 ***********************************************************************************************************
    double[][] coordinates;
    double step_size,very_small;
    
    public MinimizerException(double[][] coordinates, double step_size, double very_small){
		super();
		this.coordinates = coordinates;
		this.step_size = step_size;
		this.very_small = very_small;
    }

    public void test(Abstract_total_energy energy, Atom_list list) {
		System.out.println("Steepest_descent error:\n"+
			   "step size = "+ step_size + " < "+ very_small + "\n" + 
			   "Failed to converge. Probably a derivation problem.\n"+
			   "energy test follow:\n");
		energy.testAll(list);
		throw new RuntimeException("Steepest_descent error");
    }
*************************************************************************************************************/


