package meshi.molecularElements.allAtoms;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.excludedVol.ExcludedVolCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.tether.TetherCreator;
import meshi.energy.twoTorsions.TwoTorsionsCreator;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.ResidueCreator;
import meshi.optimizers.LBFGS;
import meshi.optimizers.LineSearchException;
import meshi.optimizers.Minimizer;
import meshi.optimizers.MinimizerException;
import meshi.optimizers.SteepestDecent;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.filters.Filter;
public class AllAtomsProtein extends Protein implements Residues , KeyWords {
    public AllAtomsProtein() {
	super();
    }
	

    public AllAtomsProtein(AtomList atoms, ResidueCreator creator) {
        super(atoms, creator);
    }

    public String relax(CommandList commands) throws MinimizerException, LineSearchException {
        System.out.println("--------- step 1 ----------------");
        CommandList minimizerCommands = commands.firstWordFilter(RELAX);
        double tolerance = minimizerCommands.secondWord(TOLERANCE).thirdWordDouble();
        System.out.println("Relaxation tolerance = "+tolerance);
        int maxIteration = minimizerCommands.secondWord(MAX_STEPS).thirdWordInt();
        System.out.println("Relaxation maxIteration = "+maxIteration);
        int reportEvery = minimizerCommands.secondWord(REPORT_EVERY).thirdWordInt();
        System.out.println("Relaxation reportEvery = "+reportEvery);
        EnergyCreator[] energyCreators2 = {new BondCreator(),
				           new AngleCreator(),
                                           new PlaneCreator(),
                                           new OutOfPlaneCreator(),
					   new TwoTorsionsCreator(),	                                      
	                                   new ExcludedVolCreator(new UpToCgamaFilter(),0.9),
					   new TetherCreator()};
        DistanceMatrix distanceMatrix = new DistanceMatrix(atoms(),5.5,0.5,4);
        TotalEnergy energy = new TotalEnergy(this, distanceMatrix, energyCreators2, commands);
        Minimizer minimizer = new SteepestDecent(energy,tolerance,maxIteration/10,reportEvery);
        minimizer.minimize();
        minimizer = new LBFGS(energy,tolerance,maxIteration,reportEvery);
        minimizer.minimize();
	System.out.println("--------- step 2 ----------------");
        EnergyCreator[] energyCreators2a = {new BondCreator(),
				           new AngleCreator(),
                                           new PlaneCreator(),
                                           new OutOfPlaneCreator(),
					    new TwoTorsionsCreator(),	                                      
					    new ExcludedVolCreator(1,1),
					    new TetherCreator()};
        distanceMatrix = new DistanceMatrix(atoms(),5.5,0.5,4);
        energy = new TotalEnergy(this, distanceMatrix, energyCreators2a, commands);
        minimizer = new SteepestDecent(energy,tolerance,maxIteration/10,reportEvery);
        minimizer.minimize();
        minimizer = new LBFGS(energy,tolerance,maxIteration,reportEvery);
        minimizer.minimize();
	System.out.println("--------- step 3 ----------------");
        EnergyCreator[] energyCreators3 = {new BondCreator(),
				           new AngleCreator(),
                                           new PlaneCreator(),
                                           new OutOfPlaneCreator(),
					   new TwoTorsionsCreator(),	                                      
					   new ExcludedVolCreator(1,1)};
        distanceMatrix = new DistanceMatrix(atoms(),5.5,0.5,4);
        energy = new TotalEnergy(this, distanceMatrix, energyCreators3, commands);
        minimizer = new SteepestDecent(energy,tolerance,maxIteration/10,reportEvery);
        minimizer.minimize();
        minimizer = new LBFGS(energy,tolerance,maxIteration,reportEvery);
        return minimizer.minimize();
    }

    //----------------------- filters ---------------------------
    private static final char[] DEZH = {'D','E','Z','H'};
    private static final char[] EZH = {'E','Z','H'};
    private static final char[] ZH = {'Z','H'};

    private class UpToCgamaFilter extends sideChainFilter {
	public UpToCgamaFilter() {
	    super(DEZH);
	}
    }

    private class UpToCdeltaFilter extends sideChainFilter {
	public UpToCdeltaFilter() {
	    super(EZH);
	}
    }

    private class UpToCepsilonFilter extends sideChainFilter {
	public UpToCepsilonFilter() {
	    super(ZH);
	}
    }

    private abstract class sideChainFilter implements Filter {
	char[] excludedSecondChars;
	int excludedSecondCharsLength;
	public sideChainFilter(char[] excludedSecondChars) {
	    this.excludedSecondChars = excludedSecondChars;
	    excludedSecondCharsLength = excludedSecondChars.length;
	}
	public boolean accept(Object obj) {
	    Distance distance = (Distance) obj;
	    String[] names = new String[2];
	    names[0] = distance.atom1().name();
	    names[1] = distance.atom2().name();
	    for (int i = 0 ;i <2; i++) {
		String name = names[i];
		if (name.length() > 1) {
		    char secondChar = name.charAt(1);
		    for (int j = 0 ; j < excludedSecondCharsLength; j++)
			if (secondChar == excludedSecondChars[j]) return false;
		}
	    }
	    return true;
	}
    }
}

