package meshi.geometry;
import java.util.Arrays;
import java.util.Iterator;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.AtomPair;
import meshi.molecularElements.AtomPairList;
import meshi.molecularElements.Protein;
import meshi.parameters.Residues;
import meshi.util.MeshiProgram;
import meshi.util.Terminator;
import meshi.util.Updateable;
import meshi.util.UpdateableException;
import meshi.util.filters.Filter;

/**
 * Where all the headache of {@link meshi.geometry.Distance distances} is handled.                   <br>
 *
 * Almost any measurable feature of a molecule is related to distances
 * between pairs of {@link meshi.molecularElements.Atom Atoms}.
 * Thus, The calculation of distances (and their inverse and derivatives)
 * are typically a computational bottleneck in structural 
 * biology applications. In all applications that we are aware of (Please,
 * enlighten us if you know better) distance calculation is done as part 
 * of the procedures that use it (say, as part of the van-der-Waals energy 
 * calculation). As a result the distance between two atoms may be calculated 
 * more then once. For example the distance between two atoms may be 
 * calculated both during angle and torsion angle energies calculations.
 * In Meshi we tried to encapsulate all distance related issues in few 
 * classes: {@link meshi.geometry.Distance Distance}, its subclasses and this one. 
 * <p> 
 * The motivation behind this class is twofold: first, to keep all the headache of 
 * cutoff distances (see below) in a single place. Second, to arrange all the                    
 *{@link meshi.geometry.Distance
 * distances                                                                                     } 
 * of a molecular system in a single data 
 * structure so that no distance is calculated more than once in a single energy 
 * evaluation. 
 *</p><p>
 * <b>
 * Distance cutoff                                                                               </b><br>
 * Calculating and storing all the distances of a molecular system requires O(n^2) time 
 * and storage, where <b>n</b> is the number of atoms. Heuristic algorithms reduces it to O(n),
 * under certain assumptions.
 * The heuristic algorithms (see references below) relays on three characteristics 
 * of energy functions and energy based simulations:
 * <ol>
 *   <li> Atoms have a characteristic minimal distance from other atoms. Thus, the distance 
 *        matrix is typically dominated by large distances.
 *   <li> Energy functions typically decay with distance and become negligible when the 
 *        atoms are far apart. <br>
 *        This implies that even if the energy function is formally defined over all distances,
 *        distances above some threshold can be considered "infinite" and having no energy contribution.
 *   <li> During most of an energy based simulation (e.g. MD, minimization, MC etc.) the atom
 *        movements are rather slow, and the distance matrix changes very little between 
 *        consecutive simulation steps.<br>
 *        This implies that most of the atom pairs that are "infinitely" distant in a given 
 *        simulation step will remain so in the next step.
 * </ol> 
 * The current implementation requires O(n^2) storage and O(n) CPU time. Future 
 * implementations are intended to be more efficient.
 * </p><p>
 *
 * The first step of the algorithm is to separate the set of distances into two groups:<br>
 * <ol>
 *   <li> bonded distances - these are distances between atoms that are not likely to be too 
 *        far apart at any time and thus, algorithms based on distance cutoffs are not applicable to them. 
 *        Typically these are atoms separated by up to three or four covalent bonds.
 *	  The number of these distances is O(n) so they do not add to the computational complexity.
 *   <li> unbonded distances - these are distances between atoms that may or may not be close in space.
 * </ol> 
 *

 * by two predefined points: rMax & rmax+buffer<br>
 * <ol>
 *    <li> 0 < distance <= rMax                                                                              <br>
 *         It is assumed that all non zero interactions occure within this distance range.                   <br>
 *         The algorithm garenties that if a pair of atoms have a distance within this range 
 *         it is included in the nonbonded list.</li>
 *    <li> rMax < distance <= rMax + buffer                                                               <br>
 *         Some of the atom pairs with a distance within this range are included in the nonbonded list.   </li>
 *    <li> rMax < distance                                                                                <br>
 *         These atom pairs are not included in the nonbonded list</li>
 *
 *  The number of atom-pairs with inter-atomic distances in the first two regions is O(n). These 
 *  atom-pairs are stored in a list called the non-bonded-list. The only O(n^2) task is to test for 
 *  each of the O(n^2) atom-pairs currently in the third region whether it moved to the first two. This 
 *  is where our assumption that changes in the inter-atomic distances are slow enters. Consider a pair of 
 *  atoms A1 and A2 that were in the coordinates C1(S) and C2(S) at some step S of the simulation
 *  such that the distance between them 
 *  D(C1(S),C2(S)) > rMax+buffer. D(C1(S+T),C2(S+T)) the distance between these atoms at some later 
 *  step S+T may be smaller or equal to rMax only if at least for one of the atoms 
 *  D(C(S),C(S+T)) > buffer/2. This condition is tested in O(n) and assuming slow atoms movements 
 *  should, for most atoms fail. Only when this condition is satisfied for one of the proteins 
 *  We should check all it's distances from other atoms again with again O(n) complexity.
 *  
 * </p><p>
 *  
 **/
public class DistanceMatrix extends MeshiProgram implements Residues, Updateable {
    public static final Terminator terminator = new Terminator();



    /*------------------------------------------ object variables --------------------------------------------*/
    public static final double DEFAULT_RMAX = 5.5;
    public static final double DEFAULT_BUFFER = 1;
    public static final int DEFAULT_BONDED_LIST_DEPTH = 5;
    private Indicator  indicatorToUpdateHB;
    private Distance infiniteDistance = new Distance();
    public Distance infiniteDistance() {return infiniteDistance;}

    /**
     * The list of all atoms in the molecular system.
     */
    protected AtomList atomList;
    /**
     * An array of all atoms in the molecular system.
     */
    protected Object[] atomArray;

    /**
     * Internal data structure.
     **/
    protected MatrixRow[] matrix;
    protected Grid grid;

    /**
     * Maximal distance for nonzero interactions.
     **/
    protected static double rMax;
    protected static double rMax2;
    protected static double edge;
    private static Double DEFAULT_EDGE(double rMax, double buffer){return rMax+2.0*buffer/3.0;}
    
    /**
     * rMax+buffer
     **/
    protected static double rMaxPlusBuffer;

    /**
     * (rMax+buffer)^2
     **/
    protected static double rMaxPlusBuffer2;

    protected double buffer;

    /**
     * (buffer/3)^2
     **/
    protected static double bufferOneThirdSqr;

    /**
     * Atom pairs with inter-atomic distances below rMax (and some of the pairs below rMax+buffer). 
     **/
    protected DistanceList nonBondedList = new DistanceList();
    
//    protected DistanceList newNonBondedList = new DistanceList();

    /**
     * List of DistanceLists
     * Every DistanceList contains distances needed for one EnergyTerm,
     * selected by its filter
     * For example:
     * - Distances of good hydrogen bonds candidate that were added in the current update opperation
     * - Applicable distances between nonBonded C-N candidate that were added in the current update opperation
     */
   protected ListOfDistanceLists energyTermsDistanceLists;
        public ListOfDistanceLists energyTermsDistanceLists(){return energyTermsDistanceLists;}

    /**
     * Atom pairs with inter-atomic distances that are always relatively small.
     **/

    protected DistanceList  bondedList;

    private int newConstant = 0;
    private boolean nonBondedFlag = true;
    // private AtomPairList addToBondedList;
    protected int bondedListDepth;
    private int numberOfUpdates = 0;
    protected boolean debug = false;
    /**
     * Enter DistanceMatrix debug mode.
     **/
    public void debugON() {debug = true;}
    /**
     * Exit DistanceMatrix debug mode.
     **/
    public void DebugOFF() {debug = false;}
    /*-------------------------------------------- constructors -----------------------------------------------*/
    protected DistanceMatrix() {}
    /**
     * A distance matrix with practicaly no distance cutoffs and with "hard" atoms.
     **/
   public DistanceMatrix(AtomList atomList) {
    this(atomList,
         DEFAULT_RMAX, //rMax
         DEFAULT_BUFFER, // buffer
         DEFAULT_EDGE(DEFAULT_RMAX,DEFAULT_BUFFER), //edge of Grid
         1// bondedListDepth
         );
    }

    public DistanceMatrix(Protein protein) {
    this(protein.atoms(), // atom list
         DEFAULT_RMAX,
         DEFAULT_BUFFER, // buffer
         DEFAULT_EDGE(DEFAULT_RMAX,DEFAULT_BUFFER), //edge of Grid
         DEFAULT_BONDED_LIST_DEPTH// bondedListDepth
         );
    }


    public DistanceMatrix(Protein protein, int bondedListDepth) {
    this(protein.atoms(), // atom list
         DEFAULT_RMAX,    //rMax
         DEFAULT_BUFFER, // buffer
         DEFAULT_EDGE(DEFAULT_RMAX,DEFAULT_BUFFER), //edge of Grid
         bondedListDepth);
    }

    public DistanceMatrix(AtomList atomList, double rMax, double buffer) {
    this(atomList,
         rMax,
         buffer,
         DEFAULT_EDGE(rMax,buffer), // edge of Grid
         DEFAULT_BONDED_LIST_DEPTH// bondedListDepth
         );
    }

    public DistanceMatrix(AtomList atomList, double rMax, double buffer, AtomPairList addToBondedList) {
    this(atomList,
         rMax, buffer,
         DEFAULT_EDGE(rMax,buffer), // edge of Grid
         DEFAULT_BONDED_LIST_DEPTH);
    }

    public DistanceMatrix(AtomList atomList, double rMax, double buffer, double edge) {
    this(atomList,
         rMax, buffer,
         edge,  // edge of Grid
         DEFAULT_BONDED_LIST_DEPTH);
    }

    public DistanceMatrix(AtomList atomList, double rMax, double buffer, int bondedListDepth) {
            this(atomList,
                  rMax,
                  buffer,
                  DEFAULT_EDGE(rMax,buffer), // edge of Grid
                  bondedListDepth);
    }
   public DistanceMatrix(AtomList atomList, double rMax, double buffer, double edge,
                         int bondedListDepth) {
    // Setting the contants used for calculating the reported distance.
    this.atomList = atomList;
    atomList.renumber();
    this.buffer = buffer;
    DistanceMatrix.rMax = rMax;
    this.edge = edge;
    //this.addToBondedList = new AtomPairList();
    this.bondedListDepth = bondedListDepth;
    reset();
    }

    private void reset() {
	terminator.reset();
	atomArray  = atomList.toArrayOfAtoms();
	Arrays.sort(atomArray);
	rMax2 = rMax*rMax;
	rMaxPlusBuffer  = rMax+buffer;
	rMaxPlusBuffer2  = rMaxPlusBuffer*rMaxPlusBuffer;
	bufferOneThirdSqr = buffer*buffer/9;
	int size = atomArray.length;
	indicatorToUpdateHB = new Indicator();
	
	// Initialize the distance matrix.
	matrix = new MatrixRow[size];
	for (int i = 0; i< size; i++)
	    matrix[i] = new MatrixRow((Atom) atomArray[i], 200, matrix);
	
	/**
	 * Every element of energyTermsDistanceLists
	 * can be specified in Creator of EnergyTerm needed it (for example HydrogenBondsCreator)
	 */
	energyTermsDistanceLists = new ListOfDistanceLists();
	
	// Generate bonded list.
	bondedList = getBondedList(atomArray, bondedListDepth, matrix);
	//bonded list.
	
	// We adjust the grid cell edge size so that the memory requirements of the
	// grid will not result in memory failure.  
	try { 
	    while ((grid = new Grid(atomArray, edge, DEFAULT_EDGE(rMax,buffer))).failToBuild()) {
		if (edge >= 100000) throw new RuntimeException("Very weird");
	       edge*=2;
	    }
	    update();
	}
	catch (UpdateableException e) {
	    throw new RuntimeException(" Cannot create grid. Apparently the atoms a spread over"+
				       " a very large volume");
	}
    }
    
    
    /*------------------------------------------ update --------------------------------------------*/
    /**
     * Updates the distance matrix. 
     **/
    public void update(int numberOfUpdates) throws UpdateableException {
    if (numberOfUpdates == this.numberOfUpdates+1) {
        this.numberOfUpdates++;
        update();
    }
    else if (numberOfUpdates != this.numberOfUpdates)
        throw new RuntimeException("Something weird with DistanceMatrix.update(int numberOfUpdates)\n"+
                       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
    }


    private void update()  throws UpdateableException{
	if (terminator.dead()) throw new RuntimeException(terminator.message());
	int length = matrix.length;
	for (int iRow = 0; iRow <length; iRow++)
	    matrix[iRow].update();
	int size = atomList.size();

        if (! grid.build()) throw new UpdateableException();
        for (int iatom = 0; iatom < size; iatom++) {
	    MatrixRow row = matrix[iatom];
	    GridCell gridCell = grid.getCell(row.atom);
	    row.addCell(gridCell);
        }
	getNonBondedList();
	//  testNonBondedList();
    }


    protected void getNonBondedList() {
	double infiniteDis = Distance.INFINITE_DISTANCE;
	int length = matrix.length;
	nonBondedList.reset();
//	newNonBondedList.reset();	
	Iterator iter = energyTermsDistanceLists.iterator();
	DistanceList distanceList;
	while ((distanceList = (DistanceList) iter.next()) != null)
	    distanceList.reset();
	
	for (int iRow = 0; iRow <length; iRow++) {
	    MatrixRow row = matrix[iRow];
        boolean doneRow = false;
        int rowLength = row.size();
        Distance[] distances = row.distances();
        for (int iDistance = 0; (! doneRow) && (iDistance < rowLength) ; iDistance++) {
        Distance distance = distances[iDistance];
        int atom2number = distance.atom2Number;
        if (atom2number >= iRow) doneRow = true;
        else
            if ((! distance.bonded) && (distance.distance < infiniteDis)){
            nonBondedList.fastAdd(distance);
            if (distance.distancePrev() == infiniteDis)   {
//                     newNonBondedList.fastAdd(distance);
                      iter = energyTermsDistanceLists.iterator();
                      while ((distanceList = (DistanceList) iter.next()) != null) 
                          distanceList.selectionAdd(distance);
                 }
             }    
        }
    }
    }

    public void testNonBondedList() {
        Distance nonBonded, distance, distanceInMatrix;
        Iterator nonBondedDistances = nonBondedList.iterator();

        for (int i = 0 ; i < atomArray.length; i++) {
            //System.out.println("Testing row "+i+" of DistanceMatrix");
            for (int j = 0; j < i; j++){
		distance = new Distance((Atom)atomArray[i], (Atom) atomArray[j]);
		if (distance.distance() < rMax) {
                    distanceInMatrix = matrix[i].binarySearch(j);
                    if ((distanceInMatrix==null) || (!distanceInMatrix.equals(distance))) {
                            System.out.println("\nError in the cell ["+i+"]["+j+"] of  DistanceMatrix");
                            System.out.println("\nDistanceMatrix Row for atom "+i+" is:");
                            matrix[i].print();
                            System.out.println("\nDistanceMatrix Row for atom "+j+" is:");
                            matrix[j].print();
                            throw new RuntimeException(" \nDistanceMatrix store the  wrong distance "+distanceInMatrix+" \ninstead the current distance "+distance);
                    }
                    else
                    if (!distanceInMatrix.bonded()) {
                        nonBonded = (Distance) nonBondedDistances.next();
                        if ((nonBonded.distance() != distance.distance()) |
                        (nonBonded.atom1() != distance.atom1()) |
                        (nonBonded.atom2() != distance.atom2())) {
                            System.out.println("Distance "+distance+" is missing from the non bonded list.");
                            System.out.println("DistanceMatrix Row for atom "+i+" is:");
                            matrix[i].print();
                            System.out.println("NonBondedList is :");
                            nonBondedList.print();
                            throw new RuntimeException("Distance "+distance+" \nis inconsistent with "+nonBonded);
                            //System.out.println("Distance "+distance+" \nis inconsistent with "+nonBonded);
                        }
//                if (nonBonded.attributes != distance.attributes)
  //                  matrix[i].print();
            //    System.out.println("Distance "+distance+" is missing from the non bonded list");
                    }
                }
            }
        }
 	}


    /*------------------------------------------ other methods --------------------------------------------*/
    public MatrixRow rowNumber(int index) { return matrix[index];}
    /**
     * Returns the non-bonded-list.
     **/
    public DistanceList nonBondedList() {
    if (nonBondedList == null)
        getNonBondedList();
    return nonBondedList;
    }
    public DistanceList bondedList() {return bondedList;}
    
//    public DistanceList newNonBondedList() {return newNonBondedList;}

     /**
     * Returns the Distance object of the parameters.
     **/
    public Distance distance(Atom atom1, Atom atom2) {
    return distance(atom1.number(), atom2.number());
    }

    public Distance distance(AtomPair atomPair) {
    return distance(atomPair.largeNumber(), atomPair.smallNumber());
    }

    /**
     * Returns the Distance object of the parameters.
     **/
    public Distance distance(int atom1Number, int atom2Number) {
    Distance out;
    out = matrix[atom1Number].binarySearch(atom2Number);
    if (out == null)
        return infiniteDistance;
    return  out;
    }

    public Distance exclusiveDistance(int atom1Number, int atom2Number) {
    return matrix[atom1Number].binarySearch(atom2Number);
    }

    public double radius() { return atomList.radius();}

    public String toString() {
    return ("DistanceMatrix:\n"+
        "\t number of atoms \t"+atomList.size()+
        "\t rMax \t"+rMax+
        "\t buffer\t"+buffer);
    }

    public String upperLeft() {
    int n = atomArray.length;
    String out = "";
    int atom1, atom2;
    Distance distance;
    if (n > 7) n = 7;
    for (int i = 0; i < n; i++) {
        atom1 = ((Atom) atomArray[i]).number();
        for (int j = 0; j < i; j++) {
        atom2 = ((Atom) atomArray[j]).number();
        distance = matrix[atom1].binarySearch(atom2);
        if (distance != null) out += distance.distance()+"\t";
        else out += "--\t";
        }
        out += "\n";
    }
    return out;
    }

    public void test() {
    int atom1, atom2;
    int n = atomArray.length;
    Distance distance;
    for (int i = 0; i < n; i++) {
        atom1 = ((Atom) atomArray[i]).number();
        for (int j = 0; j < i ; j++) {
        atom2 = ((Atom) atomArray[j]).number();
        distance =  matrix[atom1].binarySearch(atom2);
        System.out.println("DM("+atom1+","+atom2+") "+distance);
        }
    }
    Iterator distances;
    System.out.println("Bonded List:");
    distances = bondedList.iterator();
    while((distance = (Distance) distances.next()) != null)
        System.out.println((distance.atom1()+","+distance.atom2()));
    System.out.println("Non-Bonded List:");
    distances = nonBondedList.iterator();
    while((distance = (Distance) distances.next()) != null)
        System.out.println(distance.atom1()+","+distance.atom2());
    testNonBondedList();
    }
    /**
     * Returns the bonded list
     **/
    public int nonBondedListSize() {return nonBondedList.size();}

    public static DistanceList getBondedList(Object[] atomArray, int depth, MatrixRow[] matrix ) {
    AtomList bonded;
    Atom atom, bondedAtom;
    Iterator bondedAtoms;
    DistanceList out;
    int length = atomArray.length;
    AtomPairList tempList = new AtomPairList();
    for (int iatom = 0; iatom < length; iatom++) {
        atom = (Atom) atomArray[iatom];
        bonded = getBonded(atom, depth);
        bondedAtoms = bonded.iterator();
        while ((bondedAtom = (Atom) bondedAtoms.next()) != null)
        tempList.fastAdd(new AtomPair(atom, bondedAtom));
    }
    tempList.sort();
    // This instantiate the Distance objects associated with bonded atom pairs.
    // These pairs will be ignored during the non-bonded list update.
    AtomPair atomPair;
    Distance distance;
    out = new DistanceList();
    Iterator atomPairs = tempList.iterator();
    while ((atomPair = (AtomPair) atomPairs.next()) != null) {
        if (atomPair.atom1Number() > atomPair.atom2Number()) {
	    if (atomPair.atom1().frozen() & atomPair.atom2().frozen())
		distance = new FrozenDistance(atomPair.atom1(),atomPair.atom2(),
					      rMax2,rMaxPlusBuffer2);
	    else distance = new Distance(atomPair.atom1(),atomPair.atom2());
	}
        else {
	    if (atomPair.atom1().frozen() & atomPair.atom2().frozen())
		distance = new FrozenDistance(atomPair.atom2(),atomPair.atom1(),
					      rMax2,rMaxPlusBuffer2);
	    else distance = new Distance(atomPair.atom2(),atomPair.atom1());
	}
        distance.setBonded();
        //if (matrix[atomPair.largeNumber()].binarySearch(distance.atom2.number())==null)
            matrix[atomPair.largeNumber()].add(distance);
        //if (matrix[atomPair.smallNumber()].binarySearch(distance.atom2.number())==null)
            matrix[atomPair.smallNumber()].add(new DistanceMirror(distance));
        out.add(distance);
    }
    return out;
    }


    public static AtomList getBonded(Atom atom, int depth) {
    AtomList out = new AtomList();
    getBonded(atom, depth, out,atom.number());
    return out;
    }
    public static void  getBonded(Atom atom, int depth, AtomList out, int rootNumber) {
    if (depth == 0) return;
    Iterator atoms = atom.bonded().iterator();
    Atom bondedAtom;
    while ((bondedAtom = (Atom) atoms.next()) != null) {
        if ((rootNumber < bondedAtom.number()) &
        (! out.contains(bondedAtom)))
        out.fastAdd(bondedAtom);
        getBonded(bondedAtom, depth-1, out, rootNumber);
    }

    }
    public void doNotUpdateNonBondedList() {nonBondedFlag = false;}


    public static double rMax() {return rMax;}
    public static double rMax2() {return rMax2;}
    public double buffer() {return buffer;}
    public static double rMaxPlusBuffer2() {return rMaxPlusBuffer2;}
    public static double rMaxPlusBuffer() {return rMaxPlusBuffer;}
    public AtomList atomList() {return atomList;}
    public Indicator indicatorToUpdateHB() {return indicatorToUpdateHB;}

    public class Indicator{
    public Indicator(){};
    }


    private class RowIterator  implements Iterator {
    int size;
    int current;

    public RowIterator() {
        current = 0;
        size = matrix.length;
    }
    public void remove() {
        throw new RuntimeException("remove not implemented");
    }
    public Object next() {
        if (current >= size) return null;
        int temp = current;
        current++;
        return matrix[temp];
    }
    public boolean hasNext() {
        return (current < size);
    }
    }
    public Iterator rowIterator() {return new RowIterator();}
    private static class NonBondedFilter implements Filter {
        public boolean accept(Object obj) {
        Distance distance = (Distance) obj;
        boolean out = (distance.distance2 < rMax2) & (!distance.bonded());
        return out;
    }
    }



}
