package meshi.energy;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;

import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.MeshiList;
import meshi.util.Terminator;
import meshi.util.UpdateableException;
import meshi.util.Double.DoubleList;
import meshi.util.filters.Filter;
import meshi.util.formats.Fdouble;
import meshi.util.formats.Fint;
import meshi.util.formats.Format;

/**
 * A generic class for all meshi energy functions. Typically only a single such object exists, which 
 * include several energy term objects (sub classes of AbstractEnergy).
 **/
public class TotalEnergy{
	/**
	 * Number of times the energy function was evaluated. This is a simple machine-independent 
	 * measure of convergence efficiency.
	 **/
	protected int numberOfEvaluations;
	protected static int numberOfUpdates;
	protected int numberOfReports;
	protected AtomList atomList;
	protected double[][] coordinates;
	protected MeshiList energyTerms;
	protected DoubleList energyValues;
	protected double totalEnergy;
	protected Fdouble dformat = Fdouble.SHORT;
	protected Format sformat = Format.SHORT;
	protected Fint iformat = Fint.SHORT;
	protected DistanceMatrix distanceMatrix;
	private static TotalEnergy theOnlyTotalEnergy;
	public static final Terminator terminator = new Terminator();

	private static final double INFINITY = 1/0.0;

	public TotalEnergy() {}

	public TotalEnergy(AtomList atomList, 
			DistanceMatrix distanceMatrix,
			AbstractEnergy[] energyTerms) {
		this(atomList, distanceMatrix);
		for (int i = 0; i < energyTerms.length; i++)
			this.energyTerms.add(energyTerms[i]);
	}

	public TotalEnergy(Protein protein, 
			DistanceMatrix distanceMatrix,
			EnergyCreator[] energyCreators,
			CommandList commands) {
		this(protein.atoms(), distanceMatrix);
		energyTerms = new MeshiList(new IsEnergy());
		for (int i = 0; i < energyCreators.length; i++) {
			if (!energyCreators[i].weightWasSet()){
				energyCreators[i].getWeight(commands);
			}
			if (energyCreators[i].weight() != 0) {
				energyTerms.add(energyCreators[i].createEnergyTerm(protein, 
						distanceMatrix, 
						commands)); 
			}
		}
	}

	public TotalEnergy(AtomList atomList, 
			DistanceMatrix distanceMatrix) {
		this.atomList = atomList;
		coordinates = getCoordinates(atomList);

		this.distanceMatrix = distanceMatrix;
		numberOfEvaluations = 0;
		numberOfUpdates = 0;
		numberOfReports = 0;
		reset();
		terminator.reset();
		energyTerms = new MeshiList(AbstractEnergy.filter);
		resetAtomEnergies();
	}

	public TotalEnergy(AtomList atomList) {
		this(atomList,null);
	}

	public TotalEnergy(Protein protein,
			DistanceMatrix distanceMatrix,
			EnergyCreator[] energyCreators,
			CommandList commands,AtomList atoms) {
		this(atoms, distanceMatrix);
		energyTerms = new MeshiList(new IsEnergy());
		for (int i = 0; i < energyCreators.length; i++) {
			energyCreators[i].getWeight(commands);
			if (energyCreators[i].weight() != 0) {
				energyTerms.add(energyCreators[i].createEnergyTerm(protein,
						distanceMatrix,
						commands));
			}
		}
	}


	public DistanceMatrix distanceMatrix() {return distanceMatrix;}

	public void reset() {
		theOnlyTotalEnergy = this;
	}

	public void add(TotalEnergy energy) {
		atomList.add(energy.atomList);
		energyTerms.add(energy.energyTerms);
	}

	public void test() {
		// if(distanceMatrix!=null){
		// distanceMatrix.testNonBondedList();
		// }
		System.out.println("===> TotalEnergy test");
		Atom criminal = findCriminalAtom();
		if(criminal!=null){
			System.out.println("*** Criminal found:");
			System.out.println(criminal);
			AbstractEnergy criminalEnergy = findCriminalEnergyTerm(criminal);
			if(criminalEnergy != null){
				System.out.println("****************** Criminal "+criminalEnergy+" *******************");
				criminalEnergy.test(this,criminal);
			}
			else{
				for(Iterator ti = energyTerms.iterator();ti.hasNext();){
					AbstractEnergy energyTerm = (AbstractEnergy)ti.next();
					System.out.println("****************** Testing "+energyTerm+" *******************");
					energyTerm.test(this,criminal);
				}
			}
		}
		else{
			System.out.println("*** Criminal atom not founded");
			System.out.println("*** Test terminated");
		}
	}
	final static double  DX = 1e-7;
	final static String[]XYZ = new String[]{"x","y","z"};
	final static double verySmall = Math.exp(-15);

	/**
	 * Searching for a `criminal' atom
	 * @return a criminal <code>Atom</code> 
	 */
	protected Atom findCriminalAtom(){
		System.out.println("*** Criminal atom searching");

		double[][] coordinates = new double[3][];
		double x, e1 = 999999.9999999, e2 = 999999.9999999;
		double analiticalForce = 99999.99, numericalForce = -999999.99;
		double diff = 99999.9999, maxDiff = 0;
		Atom criminal=null;
		criminal_found:
			for(Iterator atoms = atomList.iterator();atoms.hasNext();){
				Atom atom = (Atom)atoms.next();
				if (! atom.frozen()) {
					coordinates[0] = atom.X();
					coordinates[1] = atom.Y();
					coordinates[2] = atom.Z();
					for(int i = 0; i< 3; i++) {
						// Whatever should be updated ( such as distance matrix torsion list etc. )
						updateDebug(); 
						x = coordinates[i][0];
						coordinates[i][1] = 0;
						e1 = evaluate();
						analiticalForce = coordinates[i][1];
						coordinates[i][0] += DX;

						updateDebug();
						e2 = evaluate();
						numericalForce = -1*(e2-e1)/DX;
						coordinates[i][0] -= DX;

						diff = Math.abs(analiticalForce - numericalForce);
						if (maxDiff < diff){// diff is maximal
							maxDiff = diff;
							System.out.println();
							System.out.println("Atom["+atom.number()+"]."+XYZ[i]+" = "+x);
							System.out.println("Analytical force = "+analiticalForce);
							System.out.println("Numerical force  = "+numericalForce);
							System.out.println("maxDiff = "+maxDiff);
							System.out.println("tolerance = "+2*maxDiff/(Math.abs(analiticalForce)+
									Math.abs(numericalForce)+verySmall));
							criminal=atom;
							// break criminal_found;
						}

						if ((e1 == AbstractEnergy.INFINITY) | (e1 == AbstractEnergy.NaN))
							System.out.println("e1 = "+e1);
						if ((e2 == AbstractEnergy.INFINITY) | (e2 == AbstractEnergy.NaN))
							System.out.println("e2 = "+e2);
						if ((analiticalForce == AbstractEnergy.INFINITY) | (analiticalForce == AbstractEnergy.NaN))
							System.out.println("analiticalForce = "+analiticalForce);
					}
					if(atom.number() != 0 && atom.number()%100 == 0)
						System.out.print(atom.number());
					else
						System.out.print(".");
				}
			}
		System.out.println();
		return criminal;
	}
	/**
	 * Searching for a `criminal' energy term , with specific atom
	 * @param atom a `criminal' <code>Atom</code>
	 * @return a criminal energy term
	 */
	protected AbstractEnergy findCriminalEnergyTerm(Atom atom){
		System.out.println("*** Criminal energy term searching");

		double[][] coordinates = new double[3][];
		double x, e1 = 999999.9999999, e2 = 999999.9999999;
		double analiticalForce = 99999.99, numericalForce = -999999.99;
		double diff = 99999.9999, maxDiff = 0;
		AbstractEnergy criminal=null;
		for(Iterator ti = energyTerms.iterator();ti.hasNext();){
			AbstractEnergy energyTerm = (AbstractEnergy)ti.next();
			coordinates[0] = atom.X();
			coordinates[1] = atom.Y();
			coordinates[2] = atom.Z();
			for(int i = 0; i< 3; i++) {
				// Whatever should be updated ( such as distance matrix torsion list etc. )
				updateDebug(); 
				x = coordinates[i][0];
				coordinates[i][1] = 0;
				e1 = energyTerm.evaluate();
				analiticalForce = coordinates[i][1];
				coordinates[i][0] += DX;

				updateDebug();
				e2 = energyTerm.evaluate();
				numericalForce = -1*(e2-e1)/DX;
				coordinates[i][0] -= DX;

				diff = Math.abs(analiticalForce - numericalForce);
				if (maxDiff < diff){// diff is maximal
					maxDiff = diff;
					System.out.println();
					System.out.println("****** test in "+energyTerm);
					System.out.println("Atom["+atom.number()+"]."+XYZ[i]+" = "+x);
					System.out.println("Analytical force = "+analiticalForce);
					System.out.println("Numerical force  = "+numericalForce);
					System.out.println("maxDiff = "+maxDiff);
					System.out.println("tolerance = "+2*maxDiff/(Math.abs(analiticalForce)+Math.abs(numericalForce)+verySmall));
					criminal= energyTerm;
				}

				if ((e1 == AbstractEnergy.INFINITY) | (e1 == AbstractEnergy.NaN))
					System.out.println("e1 = "+e1);
				if ((e2 == AbstractEnergy.INFINITY) | (e2 == AbstractEnergy.NaN))
					System.out.println("e2 = "+e2);
				if ((analiticalForce == AbstractEnergy.INFINITY) | (analiticalForce == AbstractEnergy.NaN))
					System.out.println("analiticalForce = "+analiticalForce);
			}
		}
		return criminal;
	}

	/**
	 * Finds the maximal component (in magnitude) of the gradient vecor in coordinates ( coordinates[][1] ).
	 **/
	public static double getGradMagnitude(double[][] coordinates) {
		double sumGrad = 0;
		int n = 0;
		for (int i = 0; i < coordinates.length; i++) {
			if (coordinates[i][1] != INFINITY) {
				sumGrad +=  coordinates[i][1]*coordinates[i][1];
				n++;
			}
		}
		if (n != 0) return Math.sqrt(sumGrad/n);
		return 0;
	}
	public double[][] coordinates() {return coordinates;}
	public int numberOfEvaluations() {return numberOfEvaluations;}


	protected void setCoordinates(AtomList atomList) {
		coordinates = getCoordinates(atomList);
	}

	/**
	 * Reduce the Atom coordinates of an atom list to an array.
	 * <X coordinate of atom 1><force on X coordinate of atom 1>
	 * <Y coordinate of atom 1><force on Y coordinate of atom 1>
	 * <Z coordinate of atom 1><force on Z coordinate of atom 1>
	 *                             :
	 *                             :
	 *                             :
	 *                             :
	 * <X coordinate of atom N><force on X coordinate of atom N>
	 * <Y coordinate of atom N><force on Y coordinate of atom N>
	 * <Z coordinate of atom N><force on Z coordinate of atom N>
	 **/
	private double[][] getCoordinates(AtomList atomList) {
		AtomList tempList = atomList.filter(new FreeAtom());
		int numberOfCoordinates = tempList.size()*3;
		double[][] coordinates = new double[numberOfCoordinates][];
		Iterator atomIter = tempList.iterator();
		Atom atom;
		int i = 0;
		while ((atom = (Atom) atomIter.next()) != null) {
			coordinates[i++] = atom.X();
			coordinates[i++] = atom.Y();
			coordinates[i++] = atom.Z();
		}
		return coordinates;
	}	

	protected static class FreeAtom implements Filter {
		public boolean accept(Object obj) {
			return (! ((Atom) obj).frozen());
		}
	}

	public void evaluateAtoms() {
		for (Iterator atoms = atomList.iterator(); atoms.hasNext();) {
			((Atom) atoms.next()).resetEnergy();
		}

		try {
			update();
		}
		catch (UpdateableException ex) { throw new RuntimeException(ex.toString());}

		Iterator terms = energyTerms.iterator();
		AbstractEnergy energyTerm;
		while ((energyTerm = (AbstractEnergy) terms.next()) != null) {
			energyTerm.evaluateAtoms();
		}
	}

	public double evaluate() {
		/* Currently only one instance of (a class that extends) AbstractTotalEnergy is 
	  allowed. As distanceMatrix is static */

		if (theOnlyTotalEnergy != this) 
			throw new RuntimeException("Currently only one instance of "+
					"(a class that extends) bstractTotalEnergy"+
			"is allowed. As distanceMatrix is static");  
		if (terminator.dead()) throw new RuntimeException(terminator.message());
		resetForces();
		energyValues = new DoubleList();
		double e;
		totalEnergy = 0;
		try {
			update();
		}
		catch (UpdateableException ex) {
			throw new RuntimeException("UPDATE ERROR: I added this Runtime error on 23/3/2008.");
		}
		Iterator terms = energyTerms.iterator();
		AbstractEnergy energyTerm;
		while ((energyTerm = (AbstractEnergy) terms.next()) != null) {
			e = energyTerm.evaluate();
			if ((e != 0) & (!energyTerm.isOn())) 
				throw new RuntimeException("This is weird: "+energyTerm+
						" is off but does not return zero ("+e+")");
			if ((! (e < 0) ) & (! ( e == 0)) & (!(e > 0)))System.out.println("weird energy "+energyTerm+" "+e);
			energyValues.add(new Double(e));
			totalEnergy += e;
		}
		numberOfEvaluations++;
		return totalEnergy;
	}

	/**
	 * Returns the last enrgy value that was calculated
	 **/
	public double getLastEnergy() {
		return(totalEnergy);
	}

	/**
	 * Updates all factors related to the energy function (for example the distance matrix).
	 **/
	public void update() throws UpdateableException{
		numberOfUpdates++;
		Iterator terms = energyTerms.iterator();
		AbstractEnergy energyTerm;
		UpdateableException updateableException = null;
		while ((energyTerm = (AbstractEnergy) terms.next()) != null) { 
			try { 
				energyTerm.update(numberOfUpdates);	    
			}
			catch (UpdateableException ex) {
				updateableException = ex;
			}
		} 
		if (updateableException != null) throw updateableException;
	}

	public void updateDebug() {
		if(distanceMatrix!=null){
			distanceMatrix.debugON();
		}
		try {
			update();
		}
		catch (UpdateableException ex) { throw new RuntimeException(ex.toString());}
	}



	public static double getAverageForce(double[][] coordinates) {
		double averageForce = 0;
		for (int i = 0; i < coordinates.length; i++) 
			averageForce += coordinates[i][1]*coordinates[i][1];
		averageForce = Math.sqrt(averageForce/coordinates.length);
		return averageForce;
	}

	/**
	 * Sets all forces in the coordinates to zero.
	 **/ 
	public void resetForces() {
		int length = coordinates.length;
		for (int i = 0; i < length; i++) {
			coordinates[i][1] = 0.0;
		}
	}    

	public String reportHeader() {
		String report = "";
		report += sformat.f("#")+sformat.f(" E total ")+" ";
		Iterator terms = energyTerms.iterator();
		AbstractEnergy energyTerm;
		while ((energyTerm = (AbstractEnergy) terms.next()) != null)
			report += sformat.f(energyTerm.comment())+" ";
		report += (sformat.f("maxf")+" "+sformat.f("radius")+
				" "+sformat.f("#NB")+sformat.f("#Eeval"));
		return report;
	}   

	public String report(int step) {
		for (Iterator atoms = atomList.iterator(); atoms.hasNext();) {
			Atom atom = (Atom)atoms.next();
		}


		String report = "";
		if ((numberOfReports % 10 == 0) && (step!=2)) { 
			report += reportHeader()+"\n";
		}    
		numberOfReports++;
		Iterator values = energyValues.iterator();
		Double value;
		double averageForce = getAverageForce(coordinates);
		double maxForce = getGradMagnitude(coordinates);
		report += iformat.f(step)+dformat.f(totalEnergy)+" ";
		while ((value = (Double)values.next()) != null) {
			report += dformat.f(value.doubleValue())+" ";
		}
		report += Fdouble.SHORTER.f(maxForce)+dformat.f(atomList.radius())+" ";
		if(distanceMatrix!=null){
			report += iformat.f(distanceMatrix.nonBondedListSize())+" ";
		}
		report += iformat.f(numberOfEvaluations);
		return report;
	}	

	public boolean frozenAtomsExist() {return atomList.frozenAtomsExist();}

	protected static class IsEnergy implements Filter {
		public boolean accept(Object obj) {
			return (obj instanceof AbstractEnergy);
		}
	}


	public AtomList atomList() {return atomList;}

	public DoubleList energyValues(){return energyValues;}

	/*
	 * Returns a specific energy term, according to its type
	 */
	public AbstractEnergy getEnergyTerm(AbstractEnergy ae) {
		Iterator terms = energyTerms.iterator();
		AbstractEnergy energyTerm;
		while ((energyTerm = (AbstractEnergy) terms.next()) != null) {
			if (energyTerm.getClass().equals(ae.getClass())) {
				return energyTerm;
			}
		}
		return null;
	}

	public DistanceMatrix getDistanceMatrix() {return distanceMatrix;}

	public final int numberOfUpdates() {return numberOfUpdates;}

	protected void addTerm(AbstractEnergy term){
		energyTerms.add(term);
	}

	public void resetAtomEnergies() {
		Iterator atoms = atomList.iterator();
		Atom atom;
		while ((atom = (Atom) atoms.next()) != null) {
			atom.resetEnergy();
		}
	}

	/*
	 * Returns a all instances of a specific energy term, according to its type
	 */
	public AbstractEnergy[] getEnergyTerms(AbstractEnergy ae) {
		int counter = 0;
		Iterator terms = energyTerms.iterator();
		AbstractEnergy energyTerm;
		while ((energyTerm = (AbstractEnergy) terms.next()) != null) {
			if (energyTerm.getClass().equals(ae.getClass())) {
				counter++;
			}
		}
		AbstractEnergy[] result = new AbstractEnergy[counter];
		counter = 0;
		terms = energyTerms.iterator();
		while ((energyTerm = (AbstractEnergy) terms.next()) != null) {
			if (energyTerm.getClass().equals(ae.getClass())) {
				result[counter] = energyTerm;
				counter++;
			}
		}
		return result;
	}

	public double avgEnergy() {return totalEnergy/atomList.size();}
	public MeshiList energyTerms()  {return energyTerms;}

	public void summary() {
		Iterator values = energyValues.iterator();
		Iterator terms = energyTerms.iterator();
		Double value;
		System.out.printf("\nENERGY SUMMARY\n");
		System.out.printf("%30s = %-15.8f\n","TotalEnergy",totalEnergy);
		System.out.printf("%30s = %-15.8f\n","avgEnergy",avgEnergy());

		while ((value = (Double)values.next()) != null) {
			String name = ((AbstractEnergy) terms.next()).comment();
			System.out.printf("%30s = %-15.8f\n",name,value);
		}
		SortedAtoms sortedAtoms = new SortedAtoms();
		System.out.printf("%30s = %-15.8f\n","\"cold1\" atoms energy (avg+std) = ",sortedAtoms.sum(1));
		System.out.printf("%30s = %-15.8f\n","\"cold2\" atoms energy (avg+2*std) = ",sortedAtoms.sum(2));
		System.out.printf("%30s = %-15.8f\n","average \"cold1\" atoms energy (avg+std)",
				sortedAtoms.sum(1)/sortedAtoms.lowestEnergyAtoms(1).size());
		System.out.printf("%30s = %-15.8f\n","average \"cold2\" atoms energy (avg+2*std)",
				sortedAtoms.sum(2)/sortedAtoms.lowestEnergyAtoms(2).size());
	}


	public double filteredEnergy(double stdThreshold) {
		SortedAtoms sortedAtoms = new SortedAtoms();
		return sortedAtoms.sum(stdThreshold);
	}

	public double avgFilteredEnergy(double stdThreshold) {
		SortedAtoms sortedAtoms = new SortedAtoms();
		return sortedAtoms.sum(stdThreshold)/sortedAtoms.lowestEnergyAtoms(stdThreshold).size();
	}

	public AtomList highestEnergyAtoms(double stdThreshold) {
		SortedAtoms sortedAtoms = new SortedAtoms();
		return sortedAtoms.highestEnergyAtoms(stdThreshold);
	}


	public static class EnergyComparator implements Comparator {
		public int compare(Object obj1, Object obj2) {
			Atom atom1 = (Atom) obj1;
			Atom atom2 = (Atom) obj2;
			if (atom1.energy() > atom2.energy()) return 1;
			if (atom1.energy() < atom2.energy()) return -1;
			return 0;
		}
	}


	public void off() {
		for (Iterator terms = energyTerms.iterator(); 	
		terms.hasNext();) {
			AbstractEnergy term = (AbstractEnergy) terms.next();
			term.off();
		} 
	}
	public void on() {
		for (Iterator terms = energyTerms.iterator();
		terms.hasNext();) {
			AbstractEnergy term = (AbstractEnergy) terms.next();
			term.on();
		}
	}

	private class SortedAtoms {
		private Atom[] atoms = new Atom[atomList.size()];
		double avg,std;
		public SortedAtoms() {
			for (int i = 0; i < atoms.length; i++)
				atoms[i] = (Atom) atomList().elementAt(i);
			evaluateAtoms();
			Arrays.sort(atoms,new EnergyComparator());
			double sum = 0;
			double sum2 = 0;
			double e;
			for (int i = 0; i < atoms.length; i++) {
				e = atoms[i].energy();
				sum += e;
				sum2 += e*e;
			}
			avg = sum/atoms.length;
			std = Math.sqrt(sum2/atoms.length -avg*avg);
		}
		public double avg() {return avg;}
		public double std() {return std;}


		public double sum(double stdThreshold) {
			double sum = 0;
			AtomList lowestEnergyAtoms = lowestEnergyAtoms(stdThreshold);
			Atom atom;
			for (int i = 0; i< lowestEnergyAtoms.size();i++) {
				atom = lowestEnergyAtoms.atomAt(i);
				sum+= atom.energy();
			}
			return sum;
		}

		private AtomList lowestEnergyAtoms(double stdThreshold) {
			AtomList out = new AtomList();
			int iAtom = 0;
			while ((iAtom < atoms.length) && (atoms[iAtom].energy() < avg+stdThreshold*std)) {
				out.add(atoms[iAtom]);
				iAtom++;
			}
			return out;
		}

		private AtomList highestEnergyAtoms(double stdThreshold) {
			AtomList out = new AtomList();
			int iAtom = atoms.length-1;
			while ((iAtom >= 0) &(atoms[iAtom].energy() > avg+stdThreshold*std)) {
				out.add(atoms[iAtom]);
				iAtom--;
			}
			return out;
		}
	}
}



