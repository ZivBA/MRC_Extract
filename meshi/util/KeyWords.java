package meshi.util;
public interface KeyWords {
    //---------------------------- protein -------------------------------
    final Key SEQUENCE = new Key("sequence");
    final Key SECONDARY_STRUCTURE = new Key("secondary");
    final Key MODEL_NUMBER = new Key("modelNumber");
    final Key CLASH_DISTANCE = new Key("clashDistance");
    final Key MAX_CLASHES = new Key("maxNumberOfClashes");
    final Key N_TRYS = new Key("nTrays");

    //---------------------------- minimization -------------------------------
    final Key MINIMIZE = new Key("minimize"); // First word for full minimization
    final Key RELAX = new Key("relax"); // First word for short relaxation minimization with SteepestDecent
    final Key TOLERANCE = new Key("tolerance");
    final Key MAX_STEPS = new Key("maxSteps");
    final Key REPORT_EVERY = new Key("reportEvery");
    
    //--------------------------------- energy  -----------------------------
    final Key PARAMETERS_DIRECTORY = new Key("parameters");

    //--------------------------------- energy terms -----------------------------
    final Key LENNARD_JONES_CA = new Key("LennardJonesCa");
    final Key LENNARD_JONES = new Key("LennardJones");
    final Key ANGLE_ENERGY = new Key("angleEnergy");
    final Key BOND_ENERGY = new Key("bondEnergy");
    final Key CONSENSUS_ENERGY = new Key("consensusEnergy");
    final Key PLANE_ENERGY = new Key("planeEnergy");
    final Key OUT_OFPLANE_ENERGY = new Key("outOfPlaneEnergy");
    final Key TEMPLATE_DISTANCE_CONSTRAINS = new Key("templateDistanceConstrains");
    final Key DISTANCE_CONSTRAINS_ENERGY = new Key("distanceConstrainsEnergy"); 
    final Key INFLATE_ENERGY = new Key("inflateEnergy");
    final Key VOLUME_CONSTRAIN = new Key("volumeConstrainWeight");
    final Key HYDROGEN_BONDS = new Key("hydrogenBonds");
    final Key HYDROGEN_BONDS_PAIRS = new Key("hydrogenBondsPairs");
    final Key TWO_TORSIONS_ENERGY = new Key("twoTorsionsEnergy");
    final Key FLAT_RAMACH_ENERGY = new Key("flatRamachEnergy");
    final Key PROPENSITY_TORSION_ENERGY = new Key("propensityTorsionEnergy");
    final Key ALPHA_ANGLE_ENERGY = new Key("alphaAngleEnergy");
    final Key ALPHA_TORSION_ENERGY = new Key("alphaTorsionEnergy");
    final Key SOLVATE_ENERGY = new Key("solvateEnergy");
    final Key EXCLUDED_VOL = new Key("excludedVolume");
    final Key ELECTROSTATICS = new Key("electrostatics");
    final Key DIELECTRIC_CONSTANT = new Key("dielectricConstant"); 
    //final Key EVALUATED_LOCATION_ENERGY = new Key("tetherEnergy");
    final Key TETHER_ENERGY = new Key("tetherEnergy");
    final Key CALPHA_HYDROGEN_BONDS = new Key("cAlphaHydrogenBonds");
    final Key CALPHA_HYDROGEN_BONDS_PLANE = new Key("cAlphaPlane");
    final Key HYDROGEN_BONDS_ANGLES = new Key("hydrogenBondsAngles");
    final Key HYDROGEN_BONDS_PLANE = new Key("hydrogenBondsPlane");
    //------------------------ inflate --------------------------   
    final Key RMS_TARGET = new Key("RmsTarget"); 
    //------------------------ template based distance constrains --------------------------   
    final Key INTRA_SEGMENT_FACTOR = new Key("intraSegmentFactor");
    final Key INTRA_SEGMENT_TOLERANCE = new Key("intraSegmentTolerance");
    final Key INTER_SEGMENT_FACTOR = new Key("interSegmentFactor");
    final Key INTER_SEGMENT_TOLERANCE = new Key("interSegmentTolerance");
    final Key SATURATION = new Key("saturation");
    final Key UNSATISFIED_CUTTOF = new Key("unsatisfiedCutoff");    
    final Key UP_TO_CUTOFF = new Key("upToCutoff");
    final Key DISTANCE_CONSTRAINS_MASK = new Key("constrain");

    //--------------------------------- Superimpose -------------------------------------
    final Key SUPERIMPOSE = new Key("superimpose");
    final Key REFERENCE = new Key("reference");
    final Key MODE = new Key("mode");
    final Key ALL_CA = new Key("allCa");
    

    //--------------------------------- minimization loop -------------------------------------
    final Key MINIMIZATION_LOOP = new Key("minimizationLoop");
    final Key ITERATIONS_CA = new Key("nIterationsCA");
    final Key ITERATIONS_BACKBONE = new Key("nIterationsBackbone");
    final Key ITERATIONS_ALLATOM = new Key("nIterationsAllAtoms");
    //--------------------------------- MCM -------------------------------------
    final Key MCM = new Key("MCM");
    final Key N_MCM_STEPS = new Key("numberOfMCMsteps");
    final Key INITIAL_TEMPERATURE = new Key("initialTemperature");
    final Key FINAL_TEMPERATURE = new Key("finalTemperature");

    //--------------------------------- Homology Modeling -------------------------------------
    final Key TARGET_NAME = new Key("targetName");
    final Key TEMPLATE_NAME = new Key("templatetName");
    final Key TARGET_FILE_PATH = new Key("targetFilePath");
    final Key ALINMENT_FILE_PATH = new Key("alinmentFilePath");

final Key TEMPLATE = new Key("template");

    final Key TEMPLATE_FILE_PATH = new Key("templateFilePath");
    final Key TEMPLATE_STRUCTURE = new Key("templateStructure");
    final Key TEMPLATE_DSSP = new Key("templateDssp"); 
    final Key TEMPLATE_TARGET_ALIGNMENT = new Key("templateTargetAlignment");
    final Key OUTPUT_FILE_PATH = new Key("outputFilePath");
    final Key OUTPUT_FILE_NAME = new Key("outputFileName");
    final Key SS_NAME = new Key("ssName");
    final Key LOOSEN_EDGE_LENGTH = new Key("loosenEdgeLength");
    final Key NON_FROZEN_BOND_DEPTH = new Key("nonFrozenBondDepth");
    final Key NON_FROZEN_RADIUS = new Key("nonFrozenRadius");
    final Key NUMBER_OF_MODELS = new Key("numberOfModels");

    //--------------------------------- analysis ---------------------------------
    final Key DICTIONARY_KEY = new Key("COMMENT");
    final Key MESHILOG_KEY = new Key("MESHILOG");
    final Key KEY_KEY = new Key("T1");
    final Key VALUE_KEY = new Key("V1");

    //--------------------------------- Sequencea -------------------------------------
    final Key AA_SEQUENCE = new Key("aa sequence");
    final Key SS_SEQUENCE = new Key("ss sequence");
    final Key ACCESIBILITY_SEQUENCE = new Key("accesibility sequence");
    //
    //--------------------------------- Misc -------------------------------------
    final Key ON = new Key("on");
    final Key OFF = new Key("off");
    final Key END = new Key("end");
    final Key WEIGHT = new Key("weight");
    final Key INPUT_FILE = new Key("inputFile");
    final Key CUTOFF = new Key("cutoff");
    final Key NONE = new Key("none");
    final Key USE_FAST_ARCCOS = new Key("useFastArcCos");
    final Key DRESSER_FRAGMENTS = new Key("dresserFragments");
    final Key ROTAMER_LIBRARY = new Key("rotamerLibrary");
    final Key FIX_N_TERMINAL = new Key("fixNterminal");
    final Key FIX_C_TERMINAL = new Key("fixCterminal");
}
