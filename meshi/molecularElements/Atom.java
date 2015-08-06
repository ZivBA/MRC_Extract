package meshi.molecularElements;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

import meshi.PDB.PdbLine;
import meshi.energy.TotalEnergy;
import meshi.geometry.Coordinates;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.FrozenCoordinates;
import meshi.util.Activable;
import meshi.util.Attributable;
import meshi.util.MeshiAttribute;
import meshi.util.MeshiException;
import meshi.util.MeshiList;
import meshi.util.MeshiProgram;
import meshi.util.Verbose;
import meshi.util.filters.Filter;
import meshi.util.string.StringList;
import meshi.util.string.StringParser;

/**
 * A generic atom. 
 *<br><hr>
 *<b>A possible pitfall</b> <br>
 * A {@link meshi.geometry.Coordinates Coordinates} object 
 *   ({@link meshi.molecularElements.Atom#coordinates coordinates})is associated with an atom  
 * and handle all geometric issues(position, distance, forces etc.).  
 * The "geometric" methods ({@link meshi.molecularElements.Atom#x x()},  
 *                          {@link meshi.molecularElements.Atom#y y()},   
 *                          {@link meshi.molecularElements.Atom#distanceFrom distanceFrom(Atom)}, etc.) 
 * are used to encapsulate this implementation. The {@link meshi.geometry.Coordinates Coordinates}  
 * object itself though, is publicly accessible. Removing the overhead of function calls accelerates 
 * {@link meshi.geometry.DistanceMatrix DistanceMatrix} updating which is the bottleneck of  
 * energy based calculations. <b> It may though, open the way to spectacular bugs!!!</b>. 
 *   
 * @see meshi.molecularElements.AtomList 
 **/ 										      
public class Atom  implements Comparable, Verbose, Attributable{
	//------------------------------------------------------------------------------
	//---------------------- class and object variables ----------------------------
	//------------------------------------------------------------------------------

    /**
     * The x,y, z coordinates of the atom in the Cartesian and force spaces.
     **/ 
    protected Coordinates  coordinates; 

    /**
     * The atoms name. In proteins a unique specification of it's position within the residue.
     **/
    public final String name;  

    /**
     * A unique identifier of the atom in the system. 
     **/
    public final int ID;

    private double reliability = -1;
    public void setReliability(double value) {reliability = value;}
    public double reliability() {return reliability;}
    /**
     * A (hopefully) unique identifier of the atom in the system. 
     * Among other things serves as an index into a {@link DistanceMatrix}.
     **/
    private int number;
    public final int number() {return number;}
    
    
    /**
     * Sets the atom number. Should be used with extreme caution. 
     **/
    protected void setNumber(int number) {
	this.number = number;
    }

   /**
    * The name  of the chain to which this atom belongs.
    **/
    protected String chain; 


   /**
    * Atom's Residue.
    **/
    protected Residue residue;

     /**
     * The name of the residue to which this atom belongs.
     **/
    private String residueName; 

    /**
     * The number of the residue to which this atom belongs.
     **/
    private  int residueNumber;

   /**
     * The atom type. 
     **/
    public final int type;

 
    /**
     * Atoms bonded to this atom.
     **/
    private final int ATOM_BONDED_CAPACITY = 4;
    protected AtomList bonded = new AtomList(ATOM_BONDED_CAPACITY);

    

   /**
    * Atom's Occupancy number.
    **/
    private double occupancyNumber = 1.0;

    public double occupancy() {return occupancyNumber;}
    public void setOccupancy(double occ) {occupancyNumber = occ;} 
   /**
    * Atom's Temperature factor.
    **/
    protected double temperatureFactor = -1.0;

   /**
    * Atom's Alternate location.
    **/
    protected String alternateLocation = "";

    /**
     * Number of atoms in the system.
     **/
    private static int numberOfAtoms = 0;

    /**
     * if true, this atom cannot move during minimization & MD.
     **/
    private boolean frozen = false;

    private double energy = 0.0;

    protected int EEETag = -1;

    public final boolean isCarbon;
    public final boolean isOxygen;
    public final boolean isNitrogen;
    public final boolean isSulfur;
    public final boolean isHydrogen;
    public final boolean isBackbone;
    private String pdbLine = "This atom was not generated from a PDB file";
    
    /** I added this field ('wasVisited') on 21/3/2008 for the purpose of the torsion-space minimization **/  
    public boolean wasVisited = false;
    
    /**
     * The PDB line from which this atom was built. 
     * Returns "This atom was not generated from a PDB file" no such line exists.
     *
     **/
    private String pdbLine() {return pdbLine;}
    //------------------------------------------------------------------------------
    //----------------------------- constructors -----------------------------------
    //------------------------------------------------------------------------------
    /**
     * A generic atom with unspecified coordinates.
     * Mainly for use by other more specific constructors.
     **/
    public Atom(String name,  
                Residue residue, String residueName, int residueNumber, String chain,  
                int type){
	this.name = name; 
	if (name.charAt(0) == 'C') isCarbon = true; 
	else isCarbon = false;
	if (name.charAt(0) == 'O') isOxygen = true;
	else isOxygen = false;
	if (name.charAt(0) == 'N') isNitrogen = true;
	else isNitrogen = false;
	if (name.contains("H")) isHydrogen = true;
	else isHydrogen = false;
	if (name.charAt(0) == 'S') isSulfur = true;
	else isSulfur = false;
	if (name.equals("CA") || name.equals("CB") || name.equals("N") || name.equals("C") || 
		name.equals("H") || name.equals("O") || name.equals("OXT"))
		isBackbone = true;
	else
		isBackbone = false;
	ID = numberOfAtoms;
	number = numberOfAtoms;
	numberOfAtoms++;
	this.chain = chain;
	this.residue = residue;
	this.residueName = residueName;
	this.residueNumber = residueNumber;
        this.type = type;
	coordinates = new Coordinates();
    }

    public Atom(String name, String residueName, int residueNumber, String chain, int type){
	this(name, null, residueName, residueNumber, chain, type);
    }

    public Atom(String name, String residueName, int residueNumber, int type){
	this(name, null, residueName, residueNumber, " ", type);
    }

    public Atom(String name, Residue residue, int type){
	this(name, residue, residue.name, residue.number, residue.chain, type);
    }
    
    /**
     * A generic atom based on a PDB formated line (ATOM .....).
     **/
    public Atom(PdbLine line) {
        this(line.name(),
	     line.residueName(),
	     line.residueNumber().intValue(),-1);
	coordinates.setX(line.x());
	coordinates.setY(line.y());
	coordinates.setZ(line.z());
	temperatureFactor = line.temperatureFactor();
	occupancyNumber = line.occupancyNumber();
	alternateLocation = line.alternateLocation();
	pdbLine = line.toString();
	this.chain = line.chain();
    }
 
    /**
     * A generic Atom with position specified by x, y and z.
     **/
    public Atom(double x, double y, double z,
		   String name, Residue residue, int type){
	this(name, residue, type);
	coordinates.setX(x);
	coordinates.setY(y);
	coordinates.setZ(z);
   }
 
    public Atom(double x, double y, double z,
		   String name, String residueName, int residueNumber, int type){
	this(name, residueName, residueNumber, type);
	coordinates.setX(x);
	coordinates.setY(y);
	coordinates.setZ(z);
   }
 
    /**
     * A generic Atom with random position specified by centerX, centerY, centerZ and radius.
     **/
    public Atom(double centerX, double centerY, double centerZ, double radius,
		   String name, String residueName, int residueNumber, int type){
	this(name, residueName, residueNumber, type);
	double x=radius, y=radius, z=radius;
	while (x*x + y*y + z*z > radius*radius) {
	    x = (2*MeshiProgram.randomNumberGenerator().nextDouble() - 1) * radius;
	    y = (2*MeshiProgram.randomNumberGenerator().nextDouble() - 1) * radius;
	    z = (2*MeshiProgram.randomNumberGenerator().nextDouble() - 1) * radius;
	}
	coordinates.setX(centerX+x);
	coordinates.setY(centerY+y);
	coordinates.setZ(centerZ+z);
    }
 
    public Atom(double centerX, double centerY, double centerZ, double radius,
		   String name, Residue residue, int type){
	this(name, residue.name, residue.number, type);
	double x=radius, y=radius, z=radius;
	while (x*x + y*y + z*z > radius*radius) {
	    x = (2*MeshiProgram.randomNumberGenerator().nextDouble() - 1) * radius;
	    y = (2*MeshiProgram.randomNumberGenerator().nextDouble() - 1) * radius;
	    z = (2*MeshiProgram.randomNumberGenerator().nextDouble() - 1) * radius;
	}
	coordinates.setX(centerX+x);
	coordinates.setY(centerY+y);
	coordinates.setZ(centerZ+z);
	this.residue = residue;
    }
     /**
     * A generic Atom with random position specified by centerX, centerY, centerZ, radius and min distance between atoms
     **/
  public Atom(double centerX, double centerY, double centerZ, double radius, double minDistanceSqr,
		   String name, String residueName, int residueNumber, int type, AtomList atomList){

  this(name, residueName, residueNumber, type);  
  double x, y, z, dx, dy, dz, d2;                       
  boolean newFound ;
  int nLast = atomList.size();      
  Atom atom;
   Coordinates coor;
   
   do {   newFound = true;
        do {
	    x = (2*MeshiProgram.randomNumberGenerator().nextDouble() - 1) * radius;
	    y = (2*MeshiProgram.randomNumberGenerator().nextDouble() - 1) * radius;
	    z = (2*MeshiProgram.randomNumberGenerator().nextDouble() - 1) * radius;
	} while (x*x + y*y + z*z > radius*radius) ;
            
               for (int k = 0; k<nLast; k++)
                {
                    atom = atomList.atomAt(k);
                    coor = atom.coordinates();
                    dx = coor.x()-x;
                    dy = coor.y()-y;
                    dz = coor.z()-z;
                    d2 = dx*dx+dy*dy+dz*dz;
                    if (d2 <  minDistanceSqr) { newFound = false;                                            
                                            break;}
                }                	                    
           } while (newFound == false);
  
	coordinates.setX(centerX+x);
	coordinates.setY(centerY+y);
	coordinates.setZ(centerZ+z);
    }
  
   public Atom(double centerX, double centerY, double centerZ, double radius, double minDistanceSqr,
		   String name, Residue residue, int type, AtomList atomList){
	this(name, residue.name, residue.number, type);
	  
  double x, y, z, dx, dy, dz, d2;                       
  boolean newFound ;
  int nLast = atomList.size();      
  Atom atom;
   Coordinates coor;
   
   do {   newFound = true;
        do {
	    x = (2*MeshiProgram.randomNumberGenerator().nextDouble() - 1) * radius;
	    y = (2*MeshiProgram.randomNumberGenerator().nextDouble() - 1) * radius;
	    z = (2*MeshiProgram.randomNumberGenerator().nextDouble() - 1) * radius;
	} while (x*x + y*y + z*z > radius*radius) ;
            
               for (int k = 0; k<nLast; k++)
                {
                    atom = atomList.atomAt(k);
                    coor = atom.coordinates();
                    dx = coor.x()-x;
                    dy = coor.y()-y;
                    dz = coor.z()-z;
                    d2 = dx*dx+dy*dy+dz*dz;
                    if (d2 <  minDistanceSqr) { newFound = false;                                            
                                            break;}
                }                	                    
           } while (newFound == false);
  
	coordinates.setX(centerX+x);
	coordinates.setY(centerY+y);
	coordinates.setZ(centerZ+z);
        this.residue = residue;
    }
 
    /**
     * Copy constructor.
     **/
    public Atom(Atom atom){
       this(atom.name, null , atom.residueName , atom.residueNumber , atom.chain(),
           atom.type);
       this.alternateLocation = atom.alternateLocation();
       this.occupancyNumber = atom.occupancyNumber;
       this.temperatureFactor = atom.temperatureFactor();
       this.energy = atom.energy();
       this.setXYZ(atom.x(),atom.y(),atom.z());
       chain = atom.chain;
    }

    public Atom(Residue residue, int type, Atom atom){
	this(atom.name, residue, type);
	coordinates.setX(atom.coordinates.x());
	coordinates.setY(atom.coordinates.y());
	coordinates.setZ(atom.coordinates.z());
	chain = atom.chain;
    }
    public Atom(Residue residue, int type, Atom atom, String name){
	this(name, residue, type);
	coordinates.setX(atom.coordinates.x());
	coordinates.setY(atom.coordinates.y());
	coordinates.setZ(atom.coordinates.z());
	if (atom.frozen()) freeze();
	chain = atom.chain;
    }
 
 
    //------------------------------------------------------------------------------
    //------------------------------- methods --------------------------------------
    //------------------------------------------------------------------------------

    /**
     * Returns the Atom's chain.
     * @see meshi.molecularElements.Atom#chain
     **/
    public String chain() {return chain;}
    
    /**
     * Set the atom residue.
     * @see meshi.molecularElements.Atom#residue
     **/
    public void setResidue(Residue residue) {
	if (this.residue != null) 
	    throw new MeshiException("Error while trying to set "+residue+
				  "to be the residue of "+this+"\n"+
				  "its residue is already defined as: "+this.residue);
	this.residue = residue;
    }
    
    /**
     * Returns the atom's residue.
     * @see meshi.molecularElements.Atom#residue
     **/
     public final Residue residue() {	return residue; }
     
     
    /**
     * Returns true if the protein's name is a substring of the String parameter seperated by spaces.
     * Examples: 
     * <ol><li>
     *       String s = "CA CB";            <br>
     *       Atom a;                          <br>
     *       a.name = "CB";                   <br>
     *       System.out.println(a.nameIs(s)); <br>
     *       <b> prints true </b>             <br>
     *     </li>       
     *     <li>
     *       String s = "CACB";             <br>
     *       Atom a;                          <br>
     *       a.name = "CB";                   <br>
     *       System.out.println(a.nameIs(s)); <br>
     *       <b> prints false </b>            <br>
     *      </li>
     * </ol>
     *       
     * @see meshi.molecularElements.Atom#name
     **/
    public boolean nameIs(String keys) {return nameIs(StringParser.standard(keys));}
    /**
     * Returns true if the protein's name is included in the StringList parameter.
     **/
    public boolean nameIs(StringList keys) {
	Iterator iter = keys.iterator();
	String string;
	while ((string = (String) iter.next()) != null)
	    if (string.equals(name)) return true;
	return false;
    }

    /**
     * Returns the X coordinate of the atom.
     * @see meshi.geometry.Coordinates#x
     **/
    public final double x() {return coordinates.x();}
    public final double[] X() {return coordinates.X();}
    /**
     * Returns the Y coordinate of the atom.
     * @see meshi.geometry.Coordinates#Y
     **/
    public final double y() {return coordinates.y();}
    public final double[] Y() {return coordinates.Y();}

    /**
     * Returns the Z coordinate of the atom.
     * @see meshi.geometry.Coordinates#Z
     **/
    public final double z() {return coordinates.z();}    
    public final double[] Z() {return coordinates.Z();}    

    /**
     * Sets the X coordinate of the atom.
     * @see meshi.geometry.Coordinates#setX
     **/
    public void setX(double newx) {coordinates.setX(newx);}

    /**
     * Sets the Y coordinate of the atom.
     * @see meshi.geometry.Coordinates#setY
     **/
    public void setY(double newy) {coordinates.setY(newy);}

    /**
     * Sets the Z coordinate of the atom.
     * @see meshi.geometry.Coordinates#setZ
    **/
    public void setZ(double newz) {coordinates.setZ(newz);}

    /**
     * Adds the parameter to the X coordinate of the atom.
     * @see meshi.geometry.Coordinates#addToX
    **/
    public void addToX(double addMe) {coordinates.addToX(addMe);}  

    /**
     * Adds the parameter to the Y coordinate of the atom.
     * @see meshi.geometry.Coordinates#addToY
    **/
    public void addToY(double addMe) {coordinates.addToY(addMe);}   

    /**
     * Adds the parameter to the Y coordinate of the atom.
     * @see meshi.geometry.Coordinates#addToZ
    **/
    public void addToZ(double addMe) {coordinates.addToZ(addMe);}   
  
    /**
     * Sets the X,Y,Z coordinate of the atom.
     **/
    public void setXYZ(double newx, double newy, double newz) {
	coordinates.setX(newx);
	coordinates.setY(newy);
	coordinates.setZ(newz);
    }

    /**
     * Returns the force operating on the atom in the X direction.
     * @see meshi.geometry.Coordinates#fx
     **/
    public double fx() {return coordinates.fx();} 

    /**
     * Returns the force operating on the atom in the Y direction.
     * @see meshi.geometry.Coordinates#fy
    **/
    public double fy() {return coordinates.fy();}  

    /**
     * Returns the force operating on the atom in the Z direction.
     * @see meshi.geometry.Coordinates#fz
     **/
    public double fz() {return coordinates.fz();}  

    /**
     * Sets the force operating on the atom in the X direction.
     * @see meshi.geometry.Coordinates#setFx
     **/
    public void setFx(double fx) {coordinates.setFx(fx);}  

    /**
     * Sets the force operating on the atom in the Y direction.
     * @see meshi.geometry.Coordinates#setFy
     **/
    public void setFy(double fy) {coordinates.setFy(fy);}  

    /**
     * Sets the force operating on the atom in the Z direction.
     * @see meshi.geometry.Coordinates#setFz
     **/
    public void setFz(double fz) {coordinates.setFz(fz);}  

    /**
     * Adds the parameter to the force operating on the atom in the X direction.
     * @see meshi.geometry.Coordinates#addToFx
    **/
    public void addToFx(double addMe) {coordinates.addToFx(addMe);}   

    /**
     * Adds the parameter to the force operating on the atom in the Y direction.
     * @see meshi.geometry.Coordinates#addToFy
     **/
    public void addToFy(double addMe) {coordinates.addToFy(addMe);}    

    /**
     * Adds the parameter to the force operating on the atom in the Z direction.
     * @see meshi.geometry.Coordinates#addToFz
     **/
    public void addToFz(double addMe) {coordinates.addToFz(addMe);}  

 
     /**
     * The {@link meshi.geometry.Distance Distance} between this atom and the parameter.
     **/
    public Distance distance(Atom atom) {
	return new Distance(this,atom);
    }

    /**
     * The distance between this atom and the parameter.
     **/
    public double distanceFrom(Atom atom) {
	return distance(atom).distance();
    }

    /**
     * Returns the atom as a PDB formatted String.
     **/
    public String verbose(int level) {
	String typeString;
	String frozenString;
	if (type == -1) typeString = "noType";
	else typeString = type(type);

	if (frozen()) frozenString = "frozen";
	else frozenString = "notFrozen";
	
	return toString()+" "+typeString+" "+frozenString;
    }

    public String toString() {
	return (new PdbLine(number, name, residueName,
			     chain, residueNumber,
			     x(),y(),z(), occupancyNumber, energy)
		).toString();
    }

    public String toString(int atmN) {
	return (new PdbLine(atmN, name, residueName,
			     chain, residueNumber,
			     x(),y(),z(), occupancyNumber, energy)
		).toString();
    }

    
    /** 
     * Move the atom to a random position within <b>radius</b> from (<b>centerx</b>,<b>centery</b>,<b>centerz</b>).  
     **/
    public void randomize(double radius, 
			  double centerx, double centery, double centerz){
	setX(centerx+(MeshiProgram.randomNumberGenerator().nextDouble()-0.5)*radius);
	setY(centery+(MeshiProgram.randomNumberGenerator().nextDouble()-0.5)*radius);
	setZ(centerz+(MeshiProgram.randomNumberGenerator().nextDouble()-0.5)*radius);
    }
        
    /**
     * Returns atom's comment.
     **/
    public String comment(){ return "Atom";}   

    /**
     * Returns the list of atoms bonded to this one.
     * @see meshi.molecularElements.Atom#bonded
     **/
    public AtomList bonded() {
	return bonded;
    }

    /**
     * Bonds the <b>other</b> atom  to this one.
     * Eeach atom is added to the bonded list of the other.
     **/
    public AtomPair bond(Atom other) {
	bonded.add(other);
	other.bonded.add(this);
  	return new AtomPair(other,this);
    }

    public String type() { return type(type);}
    public int iType() {return type;}
	
    public final int residueNumber() {return residueNumber;}  
    public String name() {return name;}
    public static int numberOfAtoms() {return numberOfAtoms;}

    public static void resetNumberOfAtoms() {numberOfAtoms = 0;}

    public double occupancyNumber(){ return occupancyNumber;}

    public double temperatureFactor(){ return temperatureFactor; }

    public String alternateLocation(){ return alternateLocation; }
    
    public String residueName() {return residueName;}
    
    public static int numberOfTypes() {
	return types.size();
    }
    
    public void freeze() {
	frozen = true;
	coordinates = new FrozenCoordinates(coordinates);
	TotalEnergy.terminator.kill("Some atoms froze. You must instantiate a new TotalEnergy.");
	DistanceMatrix.terminator.kill("Some atoms froze. You must instantiate a new DistanceMatrix.");
    }
    public boolean frozen() {return frozen;}
    public void defrost() {
	frozen = false;
	coordinates = new Coordinates(coordinates);
	coordinates.setFx(0);
	coordinates.setFy(0);
	coordinates.setFz(0);
	TotalEnergy.terminator.kill("Some atoms melted. You must instantiate a new TotalEnergy.");
 	DistanceMatrix.terminator.kill("Some atoms  melted. You must instantiate a new DistanceMatrix.");
   }
    public Coordinates coordinates() { return coordinates;}

    public int compareTo(Object obj) {
	Atom other = (Atom) obj;
	if (number > other.number) return 1;
	if (number < other.number) return -1;
	return 0;
    }

    /**
     * Setting atom's chain
     * @param chain a <code>String</code> value
     **/
    public void setChain(String chain) {this.chain = chain;}
    // public void shiftResidue(int shift) {residueNumber += shift;}
    // public boolean isImageAtom() {return false;}
	
    //----------------------------------------------- Atom types -------------------------------------
    /**
     * A list of atom types.
     **/
    private static MeshiList types = new MeshiList(new IsType()); //Type Objects  stored by their intValue
    private static Object[] typesArray; //Type objects  stored by associated string.

    /**
     * The list of atom types is done.
     **/
    private static boolean typesDone = false;

    private static class Type implements Comparable{
	public int type;
	public String name;
	
	public Type (String name, int type) {
	    this.name = name;
	    this.type = type;
	}
	
	public int compareTo(Object obj) {
	    return name.compareTo(((Type) obj).name);
	}
	
	public boolean equalse(Object obj) {
	    return compareTo(obj) == 0;
	}
    }
    
    public static class IsType implements Filter {
	public boolean accept(Object obj) {
	    return (obj instanceof Type);
	}
    }
    /**
     * Adds a new atom type to the types list.
     **/
    public static int addType(String name) {
	if (typesDone) throw new RuntimeException("atom type added after the list is done");
	Type type = new Type(name,types.size());
	if (types.contains(type)) throw new RuntimeException(type+" already exists in the types list");
	types.add(type);
	return type.type;
    }
    /**
     * Adds a new atom type to the types list and declare the list done.
     **/
    public static int addType(String type,String done) {
	int t;
	if (! done.equals("DONE")) 
	    throw new RuntimeException("weird parameter to addType(String type,String done): "+done);
	t = addType(type);
	typesDone = true;
	typesArray = types.toArray();
	Arrays.sort(typesArray);
	return t;
    }
    public static String type(int i) {
	return ((Type) types.elementAt(i)).name;
    }

    public static int type(String s) {
	if (s.equals("GCB")) return -1;
	int pos = Arrays.binarySearch(typesArray,new Type(s,-1));
	if (pos <= -1) {
	    throw new RuntimeException("Cannot find type "+s);
	} 
	return ((Type) typesArray[pos]).type;
    }

    /**
     * Returns the number of atom types defined.
     */
    public static int totalNumberOfAtomTypes() { return types.size(); }
    
    public void resetEnergy() {
	energy = 0;
    }

    public void addEnergy(double add) {
	energy += add;
    }

    public double energy() {return energy;}

    public void setEEETag(int tag) {EEETag = tag;}
    public int  EEETag() {return EEETag;}
    //----------------------------------- Attributes ---------------------------------------
    private HashMap attributes = new HashMap();
    public final void addAttribute(MeshiAttribute attribute){ attributes.put(attribute.key(),attribute); }
    public final MeshiAttribute getAttribute(int key){ return (MeshiAttribute)attributes.get(key); }

    // --------------------------------- for homology modeling ------------------------------
    

    public void setResidueNumber(int newNumber) {
    	if (residue==null) {
    		residueNumber = newNumber;
    	}
    	else {
    		throw new RuntimeException("\nIt is not a good idea to change the residue number of an Atom once it is linked to a residue object.\n");
    	}
    }  

    public void emptyBonded() {
            bonded = new AtomList(ATOM_BONDED_CAPACITY);
    }

    // --------------------------------- For activation/inactivation ------------------------------
    private boolean active = true;
    private Vector<Activable> activables = new Vector<Activable>();
    public void activate() {
    	active = true;
    	for (Activable activable : activables)
    		activable.updateActivity();
    }

    public void inActivate() {
    	active = false;
    	for (Activable activable : activables)
    		activable.updateActivity();    	
    }

    public void addActivable(Activable activable) {
    	activables.add(activable);
    }
    
    public boolean active() {
    	return active;
    }
   
    
}

