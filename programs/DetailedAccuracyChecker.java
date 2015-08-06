package programs;

import java.text.DecimalFormat;

import meshi.geometry.DistanceMatrix;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionList;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.KeyWords;
import meshi.util.MeshiProgram;
import meshi.util.dssp.DSSP;
import meshi.util.rotamericTools.RotamericTools;


public class DetailedAccuracyChecker extends MeshiProgram implements Residues, 
							     AtomTypes , KeyWords{ /**
									 * The implemented
									 * interfaces defines the 
									 * names of atom and residue 
									 * types. 
									 **/
    private static String solutionFileName;
    private static String modelFileName;
    private static String dsspFileName;
    private static Protein solution; 
    private static Protein model;
    private static TorsionList modelChi1,modelChi2,solutionChi1,solutionChi2; 
    private static double angleCorrect;
    private static DistanceMatrix distanceMatrix;
    private static DSSP dssp;
    private static double TH_relACC = 0.2;
    private static String detailRes = "";
    private static DecimalFormat fmt1 = new DecimalFormat("0.#");    
    private static DecimalFormat fmt2 = new DecimalFormat("0.##");    
    
    
    public static void main(String[] args) {
    int[] nRes = new int[20];
    int[] chi1Correct = new int[20];
    int[] chi1chi2Correct = new int[20];
    int[] nResCore = new int[20];
    int[] chi1CoreCorrect = new int[20];
    int[] chi1chi2CoreCorrect = new int[20];
    int nResChi1Tot,chi1CorrectTot,nResChi2Tot,chi2CorrectTot,nResChi1CoreTot,chi1CorrectCoreTot,nResChi2CoreTot,chi2CorrectCoreTot;
    int resNum,resType;
    
    double[] chi1Diff = new double[20];
    double[] chi2Diff = new double[20];
    double[] chi1DiffCore = new double[20];
    double[] chi2DiffCore = new double[20];
    double[] rmsd = new double[20];
    double[] rmsdCore = new double[20];
    double chi1DiffTot,chi2DiffTot,chi1DiffCoreTot,chi2DiffCoreTot,rmsdTot,rmsdCoreTot;
    
    for (int c=0 ; c<20 ; c++) {
    	nRes[c] = nResCore[c] = chi1Correct[c] = chi1CoreCorrect[c] = 
    	chi1chi2CoreCorrect[c] = chi1chi2Correct[c] = 0;
    	chi1Diff[c] = chi2Diff[c] = chi1DiffCore[c] = chi2DiffCore[c] = rmsd[c] = rmsdCore[c] = 0.0;
    }
  	nResChi1Tot = chi1CorrectTot = nResChi2Tot = chi2CorrectTot = nResChi1CoreTot = chi1CorrectCoreTot = nResChi2CoreTot = chi2CorrectCoreTot = 0;
  	chi1DiffTot = chi2DiffTot = chi1DiffCoreTot = chi2DiffCoreTot = rmsdTot = rmsdCoreTot = 0.0;
	
    Torsion tor1s,tor1m,tor2s,tor2m;
    double diff,diffalt,rms,relACC;

	init(args); 
	try {
	solution = new Protein((new AtomList(solutionFileName)).filter(new AtomList.NoAlternativeLocationFilter()),
		                    new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
	model = new Protein((new AtomList(modelFileName)),new ResidueExtendedAtoms(DO_NOT_ADD_ATOMS));
	dssp = new DSSP(dsspFileName);
    }
    catch (Exception e) {
    	System.out.print("Can not find file(s)\n"+e);
    	System.exit(0);
    }
  	
  	// Building the torsion lists.
  	distanceMatrix = new DistanceMatrix(model.atoms(), 6.0, 2.0,4); 
  	try {
  		distanceMatrix.update(1); 
  	}
    catch (Exception e) {
    	throw new RuntimeException(e);
    }
    modelChi1 = (TorsionList) TorsionList.createTorsionList(model,distanceMatrix).filter(new TorsionList.FilterChi1() , new TorsionList());
    modelChi2 = (TorsionList) TorsionList.createTorsionList(model,distanceMatrix).filter(new TorsionList.FilterChi2() , new TorsionList());
  	distanceMatrix = new DistanceMatrix(solution.atoms(), 6.0, 2.0,4); 
  	try {
  		distanceMatrix.update(1); 
  	}
    catch (Exception e) {
    	throw new RuntimeException(e);
    }
    solutionChi1 = (TorsionList) TorsionList.createTorsionList(solution,distanceMatrix).filter(new TorsionList.FilterChi1() , new TorsionList());
    solutionChi2 = (TorsionList) TorsionList.createTorsionList(solution,distanceMatrix).filter(new TorsionList.FilterChi2() , new TorsionList());
    
    // Actual analisys of the torsion lists
    for (int c=0 ; c<solutionChi1.size() ; c++) {
    	tor1s = solutionChi1.torsionAt(c);
    	resNum = tor1s.getTorsionResNum();
    	resType = Residue.type(tor1s.getTorsionResName());
    	tor1m = getTor(modelChi1,resNum); 
    	tor2s = getTor(solutionChi2,resNum); 
    	tor2m = getTor(modelChi2,resNum); 
    	rms = RotamericTools.calcRMS(solution.residue(resNum),model.residue(resNum));
    	relACC = dssp.relACCofRes(resNum,' ');
    	if ((rms>=0) && (relACC>=0) && (tor1m!=null) && 
    	   ((resType==1) || (resType==15) || (resType==16) || (resType==17) ||
    	    ((tor2m!=null) && (tor2s!=null)))) {
    	    	// Thing you do anyway...
    	    	nRes[resType]++;
    	    	nResChi1Tot++;
    	        if ((tor2m != null) && (tor2s != null)) { 
    			    nResChi2Tot++;
   			    	if (relACC<TH_relACC) {
   			    		nResChi2CoreTot++;
   			    	}
   			    }
    	    	if (relACC<TH_relACC) {
    	    		nResCore[resType]++;
    	    		nResChi1CoreTot++;
    	    		rmsdCore[resType] += rms;
    	    		rmsdCoreTot += rms;
    	    	}
    	    	rmsd[resType] += rms;
    	    	rmsdTot += rms;
    	    	diff = Math.abs(tor1m.torsion()-tor1s.torsion());
    			if (Math.abs(diff-2*Math.PI) < diff)
    			   diff = Math.abs(diff-2*Math.PI);
    			// Thing you do when chi1 is correct.
    			if (diff < angleCorrect) {
    				chi1Correct[resType]++;
    				chi1CorrectTot++;
    				if (relACC<TH_relACC) {
    					chi1CoreCorrect[resType]++;
    					chi1CorrectCoreTot++;
    					chi1DiffCore[resType] += diff;
    					chi1DiffCoreTot += diff;
    				}
    				chi1Diff[resType] += diff;
    				chi1DiffTot += diff;
    				// Checking chi2...    				
    			    if ((tor2m != null) && (tor2s != null)) {
    			    	diff = Math.abs(tor2m.torsion()-tor2s.torsion());
    			    	if (tor2m.torsion()<0)
    			    	   diffalt = Math.abs(tor2m.torsion()+Math.PI-tor2s.torsion());
    			    	else
    			    	   diffalt = Math.abs(tor2m.torsion()-Math.PI-tor2s.torsion());
    			    	if (Math.abs(diff-2*Math.PI) < diff)
    			    	   diff = Math.abs(diff-2*Math.PI);
    			    	if (Math.abs(diffalt-2*Math.PI) < diffalt)
    			    	   diffalt = Math.abs(diffalt-2*Math.PI);
    			        if ((resType==2) || (resType==4) || (resType==19))
                           diff = Math.min(diff,diffalt);
    			    	if (diff < angleCorrect) {
    			    	   chi1chi2Correct[resType]++;
    			    	   chi2CorrectTot++;
    			    	   chi2Diff[resType]+=diff;
    			    	   chi2DiffTot+=diff;
    			    	   if (relACC<TH_relACC) {
    			    	   	   chi1chi2CoreCorrect[resType]++;
    			    	   	   chi2CorrectCoreTot++;
    			    	   	   chi2DiffCore[resType]+=diff;
    			    	   	   chi2DiffCoreTot+=diff;    			    	   	   
    			    	   }
    			    	}
    			    	else if (resType==Residue.type(detailRes)) 
    			    		System.out.println("BAD " + detailRes + "(2): " + tor1s.getTorsionResNum());
    			    }
    	    	}
    	    	else if (resType==Residue.type(detailRes)) 
    	    		System.out.println("BAD " + detailRes + "(1): " + tor1s.getTorsionResNum());  
    	    }
    }
    for (int c=0 ; c<20 ; c++) 
    	System.out.println(Residue.nameThreeLetters(c) + " " + nRes[c]
    	 + " " +  chi1Correct[c] + " " +  chi1chi2Correct[c] + " " +  
    	chi1Diff[c] + " " + chi2Diff[c] + " " + rmsd[c] + "      " + nResCore[c]
    	 + " " +  chi1CoreCorrect[c] + " " +  chi1chi2CoreCorrect[c] + " " +  
    	 rmsdCore[c] + " " + chi1DiffCore[c] + " " + chi2DiffCore[c]);
    	System.out.println("ALL  " + " " + nResChi1Tot
    	 + " " +  fmt1.format(chi1CorrectTot*100.0/nResChi1Tot)
    	 + " " +  fmt1.format(chi2CorrectTot*100.0/nResChi2Tot) + " " +  
    	fmt1.format(180.0/Math.PI*chi1DiffTot/chi1CorrectTot) + " " + 
    	fmt1.format(180.0/Math.PI*chi2DiffTot/chi2CorrectTot) + " " + 
    	fmt2.format(Math.sqrt(rmsdTot/nResChi1Tot)) + "      " + 
    	nResChi1CoreTot
    	 + " " +  fmt1.format(chi1CorrectCoreTot*100.0/nResChi1CoreTot)
    	 + " " +  fmt1.format(chi2CorrectCoreTot*100.0/nResChi2CoreTot) + " " +  
    	 fmt1.format(180.0/Math.PI*chi1DiffCoreTot/chi1CorrectCoreTot) + " " + 
    	 fmt1.format(180.0/Math.PI*chi2DiffCoreTot/chi2CorrectCoreTot) + " " + 
    	 fmt2.format(Math.sqrt(rmsdCoreTot/nResChi1CoreTot)));
    }    
    
    
    /** ================================= init =========================================
     *
     *A static function for parsing of the command line arguments and assigning the 
     *variables commandsFileName, modelFileName and randomNumberSeed with the right inputs. Note that this
     *static method is using parsing functions such as getOrderedArguments that are defined in MeshiProgram
     *that MinimizeProtein inherits.
     **/
     
    protected static void init(String[] args) {
 
	/**** NOTE *** the next two lines. Because of a BUG in the Java VM, the 
	 * interfaces "Residues" and "AtomTypes" are not loaded automatically when MinimizeProtein initialize. 
	 * For this purpose these two lines are crucial wherever these two interfaces are implemented. The user might 
	 * rightfully feel that these two lines are "black magic" programming, but happily to our knowledge this is 
	 * the only bizarre phenomenon we are aware of in meshi.
	 **/
	int zvl = ALA; // force the reading of "meshi.parameters.Residues"
	zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"


	String errorMessage = ("\n                  ******************\n"+
			       "Usage java -Xmx300m AccuracyChecker <solution pdb> <model pdb> <correct angle>\n"+
			       "                    ******************\n");
			      
	if (getFlag("-debug",args)) tableSet("debug",new Boolean(true));

	solutionFileName = getOrderedArgument(args);
	if (solutionFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# solution file name is "+solutionFileName);
	
	modelFileName = getOrderedArgument(args);
	if (modelFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# model file name is "+modelFileName);

	dsspFileName = getOrderedArgument(args);
	if (dsspFileName == null) throw new RuntimeException(errorMessage);
	System.out.println("# dssp file name is "+dsspFileName);

	String angleCorrectString = getOrderedArgument(args);
	if (angleCorrectString== null) throw new RuntimeException(errorMessage);
	angleCorrect = (new Double(angleCorrectString)).doubleValue();
	System.out.println("# Correct angles are under "+angleCorrect);
	angleCorrect = angleCorrect*Math.PI/180.0;

	String accessString = getOrderedArgument(args);
	if (accessString!= null) {
		TH_relACC = (new Double(accessString)).doubleValue();
		System.out.println("# Relative accessability "+TH_relACC);
	}

	String detailString = getOrderedArgument(args);
	if (detailString!= null) {
		detailRes = detailString.trim();
	}

	initRandom(0);
    }

    private static Torsion getTor(TorsionList list, int resNum) {
    	for (int cc=0 ; cc<list.size() ; cc++) 
    		if (list.torsionAt(cc).getTorsionResNum()==resNum)
    		   return list.torsionAt(cc);
    	return null;
    }




}
