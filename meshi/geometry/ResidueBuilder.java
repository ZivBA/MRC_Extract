package meshi.geometry;
import java.util.Iterator;

import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Residue;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;


public class ResidueBuilder implements Residues,AtomTypes {
    
    private static double[] xyza = new double[4];

    public ResidueBuilder() {
	    int zvl = ALA; // force the reading of "meshi.parameters.Residues"
	    zvl = ACA;// force the reading of "meshi.parameters.AtomTypes"
    }

    public static void build(Residue res, int type, double[] chi1to4) {
        switch (type){
        	
        case 0:  buildALA( res) ;//A
                 break;
        case 1:  buildCYS( res, chi1to4) ;//C
                 break;
        case 2:  buildASP( res, chi1to4) ;//D
                 break;
        case 3:  buildGLU( res, chi1to4) ;//E
                 break;
        case 4:  buildPHE( res, chi1to4) ;//F
                 break;
        case 5:  buildGLY( res) ;//G
                 break;
        case 6:  buildHIS( res, chi1to4) ;//H
                 break;
        case 7:  buildILE( res, chi1to4) ;//I
                 break;
        case 8:  buildLYS( res, chi1to4) ;//K
                 break;
        case 9:  buildLEU( res, chi1to4);//L
                 break;
        case 10:  buildMET( res, chi1to4) ;//M
                 break;
        case 11:  buildASN( res, chi1to4) ;//N
                 break;
        case 12:  buildPRO( res, chi1to4) ;//P
                 break;
        case 13:  buildGLN( res, chi1to4) ;//Q
                 break;
        case 14:  buildARG( res, chi1to4) ;//R
                 break;
        case 15:  buildSER( res, chi1to4) ;//S
                 break;
         case 16:  buildTHR( res, chi1to4) ;//T
                 break;
         case 17:  buildVAL( res, chi1to4) ;//V
                 break;
         case 18:  buildTRP( res, chi1to4) ;//W
                 break;
         case 19:  buildTYR( res, chi1to4) ;//Y
                 break;
        default: throw new RuntimeException("wrong residue type: "+type);
        }
    }

    /*
     * Putting the coordinates (in xyza) of the forth atom in a torsion, given the other three.
     */
    public static void  getAtom_xyza(double[] xyza ,double bond,double angle,double chi,Atom A3,Atom A2,Atom A1){
    	if (A3==null) 
    		throw new RuntimeException("Missing atom 3 in the torsion construction");
    	if (A2==null) 
    		throw new RuntimeException("Missing atom 2 in the torsion construction");
    	if (A1==null)
    		throw new RuntimeException("Missing atom 1 in the torsion construction");
        /*xyza should be of type double[4]*/
        for (int i=0;i<4;i++)
            xyza[i] =0;
        double moveTo[] = new double[4];
        moveTo[0] = bond*Math.sin(angle)*Math.sin(chi);
        moveTo[1] =bond*Math.sin(angle)*Math.cos(chi);
        moveTo[2] =  bond*Math. cos(angle);
        moveTo[3] = 1;

        double M[][] = new double[4][4] ;
        double invM[][] = new double[4][4] ;
        ViewAt.transformToOrigin(M,invM,A3, A2,A1);
        // moveTo*invM
        for (int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                xyza[i] += moveTo[j] * invM[j][i];
            }
        }
    }


    public static void buildBackbone(Residue originalRes,double phi,double psi) {
    	// Build O
        double bond=1.233;     
        double angle=(Math.PI/180)*120.54; 

        getAtom_xyza(xyza,  bond, angle, psi+Math.PI,
        getAtom(originalRes.atoms() , "C"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "O").setXYZ(xyza[0],  xyza[1], xyza[2]);    	
        
        // Build H
        if (getAtom(originalRes.atoms() , "H")==null)
        	return;
        bond=1.0;     
        angle=(Math.PI/180)*119.0; 

        getAtom_xyza(xyza,  bond, angle, phi+Math.PI,
        getAtom(originalRes.atoms() , "N"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "C"));
        
        getAtom(originalRes.atoms() , "H").setXYZ(xyza[0],  xyza[1], xyza[2]);    	
        
    }
    
    public static void buildBackboneFromNterm(Residue prevRes, Residue buildRes, double prevPsi, double omg,  double phi,double psi) {
       	// Build N
        double bond=1.33;     
        double angle=(Math.PI/180)*116.8; 
        getAtom_xyza(xyza,  bond, angle, prevPsi,
        getAtom(prevRes.atoms() , "C"), getAtom(prevRes.atoms() , "CA"), getAtom(prevRes.atoms() , "N"));
        getAtom(buildRes.atoms() , "N").setXYZ(xyza[0],  xyza[1], xyza[2]);    	
    	
       	// Build CA
        bond=1.46;     
        angle=(Math.PI/180)*121.3; 
        getAtom_xyza(xyza,  bond, angle, omg,
        getAtom(buildRes.atoms() , "N"), getAtom(prevRes.atoms() , "C"), getAtom(prevRes.atoms() , "CA"));
        getAtom(buildRes.atoms() , "CA").setXYZ(xyza[0],  xyza[1], xyza[2]);    	
    	
        // Build H
        if (getAtom(buildRes.atoms() , "H")!=null) {
        	bond=1.0;     
        	angle=(Math.PI/180)*119.3; 
        	getAtom_xyza(xyza,  bond, angle, Math.PI,
        			getAtom(buildRes.atoms() , "N"), getAtom(prevRes.atoms() , "C"), getAtom(prevRes.atoms() , "O"));
        	getAtom(buildRes.atoms() , "H").setXYZ(xyza[0],  xyza[1], xyza[2]);    	
        }
    	
       	// Build C
        bond=1.53;     
        angle=(Math.PI/180)*111.0; 
        getAtom_xyza(xyza,  bond, angle, phi,
        getAtom(buildRes.atoms() , "CA"), getAtom(buildRes.atoms() , "N"), getAtom(prevRes.atoms() , "C"));
        getAtom(buildRes.atoms() , "C").setXYZ(xyza[0],  xyza[1], xyza[2]);    	
        
       	// Build O
        bond=1.23;     
        angle=(Math.PI/180)*120.5; 
        getAtom_xyza(xyza,  bond, angle, psi + Math.PI,
        getAtom(buildRes.atoms() , "C"), getAtom(buildRes.atoms() , "CA"), getAtom(buildRes.atoms() , "N"));
        getAtom(buildRes.atoms() , "O").setXYZ(xyza[0],  xyza[1], xyza[2]);    	
    }
    
    
    public static void buildBackboneFromCterm(Residue postRes, Residue buildRes, double phi, double psi,  double postOmg, double postPhi) {
       	// Build C
        double bond=1.33;     
        double angle=(Math.PI/180)*121.3; 
        getAtom_xyza(xyza,  bond, angle, postPhi,
        getAtom(postRes.atoms() , "N"), getAtom(postRes.atoms() , "CA"), getAtom(postRes.atoms() , "C"));
        getAtom(buildRes.atoms() , "C").setXYZ(xyza[0],  xyza[1], xyza[2]);    	
    	
       	// Build CA
        bond=1.53;     
        angle=(Math.PI/180)*116.8; 
        getAtom_xyza(xyza,  bond, angle, postOmg,
        getAtom(buildRes.atoms() , "C"), getAtom(postRes.atoms() , "N"), getAtom(postRes.atoms() , "CA"));
        getAtom(buildRes.atoms() , "CA").setXYZ(xyza[0],  xyza[1], xyza[2]);    	

       	// Build O
        bond=1.23;     
        angle=(Math.PI/180)*122.6; 
        getAtom_xyza(xyza,  bond, angle, postOmg+Math.PI,
        getAtom(buildRes.atoms() , "C"), getAtom(postRes.atoms() , "N"), getAtom(postRes.atoms() , "CA"));
        getAtom(buildRes.atoms() , "O").setXYZ(xyza[0],  xyza[1], xyza[2]);    	

        // Build N
        bond=1.46;     
        angle=(Math.PI/180)*111.0; 
        getAtom_xyza(xyza,  bond, angle, psi,
        getAtom(buildRes.atoms() , "CA"), getAtom(buildRes.atoms() , "C"), getAtom(postRes.atoms() , "N"));
        getAtom(buildRes.atoms() , "N").setXYZ(xyza[0],  xyza[1], xyza[2]);    	
    	
        // Build H
        if (getAtom(buildRes.atoms() , "H")!=null) {
        	bond=1.0;     
        	angle=(Math.PI/180)*119.0; 
        	getAtom_xyza(xyza,  bond, angle, phi + Math.PI,
        			getAtom(buildRes.atoms() , "N"), getAtom(buildRes.atoms() , "CA"), getAtom(buildRes.atoms() , "C"));
        	getAtom(buildRes.atoms() , "H").setXYZ(xyza[0],  xyza[1], xyza[2]);    	
        }        
    }

    protected static void buildALA(Residue originalRes) {
    	// Build CB
        double bond=1.54;     
        double angle=(Math.PI/180)*109.8; 

        getAtom_xyza(xyza,  bond, angle, 122.2*(Math.PI/180),
        getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "C"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CB").setXYZ(xyza[0],  xyza[1], xyza[2]);    	
    }


    protected static void buildCYS(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
        //create SG
        double bond=1.808;     /*        CH2E-SG1E */
        double angle=(Math.PI/180)*114.4; /*  CH1E-CH2E-SG1E */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "SG").setXYZ(xyza[0],  xyza[1], xyza[2]);
    }

    protected static void buildASP(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
        /* Atom ASP CG C */
        double bond=1.516;     /*       CH2E-C */
        double angle=(Math.PI/180)*112.6; /*  CH1E-CH2E-C */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom ASP OD1 OC */
        bond= 1.249;     /*       C-OC */
        angle=(Math.PI/180)*118.4;/*  CH2E-C-OC */

        getAtom_xyza(xyza,  bond, angle, chi[1],
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "OD1").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom ASP OD2 OC */
        bond= 1.249;     /*       C-OC */
        angle=(Math.PI/180)*118.4;/*  CH2E-C-OC */

        getAtom_xyza(xyza,  bond, angle, chi[1]-Math.PI,
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "OD2").setXYZ(xyza[0],  xyza[1], xyza[2]);
    }

    protected static void buildGLU(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
        /* Atom GLU CG CH2E */
        double bond=1.520;     /*       CH2E-CH2E */
        double angle=(Math.PI/180)*114.1; /*  CH1E-CH2E-CH2E */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom GLU CD C */
        bond= 1.516;     /*       CH2E-C */
        angle=(Math.PI/180)*112.6; /*  CH2E-CH2E-C */

        getAtom_xyza(xyza,  bond, angle, chi[1],
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "CD").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom GLU OE1 OC */
        bond= 1.249;     /*       C-OC */
        angle=(Math.PI/180)*118.4; /*  CH2E-C-OC */

        getAtom_xyza(xyza,  bond, angle, chi[2],
        getAtom(originalRes.atoms() , "CD"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
        
        getAtom(originalRes.atoms() , "OE1").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom GLU OE2 OC */
        bond= 1.249;     /*       C-OC */
        angle=(Math.PI/180)*118.4; /*  CH2E-C-OC */

        getAtom_xyza(xyza,  bond, angle, chi[2]+Math.PI,
        getAtom(originalRes.atoms() , "CD"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
        
        getAtom(originalRes.atoms() , "OE2").setXYZ(xyza[0],  xyza[1], xyza[2]);
    }

    protected static void buildPHE(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
        /* Atom PHE CG CF */
        double bond=1.502;     /*       CH2E-CF */
        double angle=(Math.PI/180)*113.8; /*  CH1E-CH2E-CF */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom PHE CD1 CR1E */
        bond=1.384;     /*       CF-CR1E */
        angle=(Math.PI/180)*120.7; /*  CH2E-CF-CR1E */

        getAtom_xyza(xyza,  bond, angle, chi[1],
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "CD1").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom PHE CD2 CR1E */
        bond=1.384;     /*       CF-CR1E */
        angle=(Math.PI/180)*120.7; /*  CH2E-CF-CR1E */

        getAtom_xyza(xyza,  bond, angle, chi[1]-Math.PI,
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "CD2").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom PHE CE1 CR1E */
        bond=1.382;     /*     CR1E-CR1E */
        angle=(Math.PI/180)*120.7; /*  CF-CR1E-CR1E */

        getAtom_xyza(xyza,  bond, angle, Math.PI,
        getAtom(originalRes.atoms() , "CD1"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
        
        getAtom(originalRes.atoms() , "CE1").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom PHE CE2 CR1E */
        bond=1.382;     /*     CR1E-CR1E */
        angle=(Math.PI/180)*120.7; /*  CF-CR1E-CR1E */

        getAtom_xyza(xyza,  bond, angle, Math.PI,
        getAtom(originalRes.atoms() , "CD2"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
        
        getAtom(originalRes.atoms() , "CE2").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom PHE CZ CR1E */
        bond=1.382;     /*       CR1E-CR1E */
        angle=(Math.PI/180)*120.0; /*  CR1E-CR1E-CR1E */

        getAtom_xyza(xyza,  bond, angle, 0,
        getAtom(originalRes.atoms() , "CE2"), getAtom(originalRes.atoms() , "CD2"), getAtom(originalRes.atoms() , "CG"));
        
        getAtom(originalRes.atoms() , "CZ").setXYZ(xyza[0],  xyza[1], xyza[2]);
    }

    protected static void buildGLY(Residue originalRes) {    }


    protected static void buildHIS(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
        /* Atom HIS CG C5 */
       double bond=1.497;     /*       CH2E-C5 */
       double angle=(Math.PI/180)*113.8; /*  CH1E-CH2E-C5 */

       getAtom_xyza(xyza,  bond, angle, chi[0],
       getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
       
       getAtom(originalRes.atoms() , "CG").setXYZ(xyza[0],  xyza[1], xyza[2]);


       /* Atom HIS ND1 NH1 */
       bond=1.378;     /*       C5-NH1 */
       angle=(Math.PI/180)*122.7; /*  CH2E-C5-NH1 */

       getAtom_xyza(xyza,  bond, angle, chi[1],
       getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
       
       getAtom(originalRes.atoms() , "ND1").setXYZ(xyza[0],  xyza[1], xyza[2]);

       /* Atom HIS HD1 NH1 */
       bond=1.0;     /*       C5-NH1 */
       angle=(Math.PI/180)*125; /*  CH2E-C5-NH1 */

       getAtom_xyza(xyza,  bond, angle, 0,
       getAtom(originalRes.atoms() , "ND1"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
       
       getAtom(originalRes.atoms() , "HD1").setXYZ(xyza[0],  xyza[1], xyza[2]);

       /* Atom HIS CE1 CRH */
       bond=1.345;     /*     NH1-CRH */
       angle=(Math.PI/180)*109.0; /*  C5-NH1-CRH */

       getAtom_xyza(xyza,  bond, angle, Math.PI,
       getAtom(originalRes.atoms() , "ND1"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
       
       getAtom(originalRes.atoms() , "CE1").setXYZ(xyza[0],  xyza[1], xyza[2]);

       /* Atom HIS NE2 NR */
       bond=1.319;     /*      CRH-NR */
       angle=(Math.PI/180)*111.7; /*  NH1-CRH-NR */

       getAtom_xyza(xyza,  bond, angle, 0,
       getAtom(originalRes.atoms() , "CE1"), getAtom(originalRes.atoms() , "ND1"), getAtom(originalRes.atoms() , "CG"));
       
       getAtom(originalRes.atoms() , "NE2").setXYZ(xyza[0],  xyza[1], xyza[2]);

       /* Atom HIS HE2 */
       bond=1.0;     /*       C5-NH1 */
       angle=(Math.PI/180)*125.57; /*  CH2E-C5-NH1 */

       getAtom_xyza(xyza,  bond, angle, Math.PI,
       getAtom(originalRes.atoms() , "NE2"), getAtom(originalRes.atoms() , "CE1"), getAtom(originalRes.atoms() , "ND1"));
       
       getAtom(originalRes.atoms() , "HE2").setXYZ(xyza[0],  xyza[1], xyza[2]);

       /* Atom HIS CD2 CR1E */
       bond=1.356;     /*      NR-CR1E */
       angle=(Math.PI/180)*130.69; /*  CRH-NR-CR1E */

       getAtom_xyza(xyza,  bond, angle, chi[1]-Math.PI,
       getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
       
       getAtom(originalRes.atoms() , "CD2").setXYZ(xyza[0],  xyza[1], xyza[2]);
    }

    protected static void buildILE(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
    /* Atom ILE CG1 CH2E */
   double bond=1.530;     /*       CH1E-CH2E */
   double angle=(Math.PI/180)*110.4; /*  CH1E-CH1E-CH2E */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG1").setXYZ(xyza[0],  xyza[1], xyza[2]);

    /* Atom ILE CG2 CH3E */
    bond=1.521;     /*       CH1E-CH3E */
    angle=(Math.PI/180)*110.5; /*  CH1E-CH1E-CH3E */

        getAtom_xyza(xyza,  bond, angle, chi[0]- (Math.PI*2.0/3.0),
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG2").setXYZ(xyza[0],  xyza[1], xyza[2]);

    /* Atom ILE CD CH3E */
    bond=1.513;     /*       CH2E-CH3E */
    angle=(Math.PI/180)*113.8; /*  CH1E-CH2E-CH3E */

        getAtom_xyza(xyza,  bond, angle, chi[1],
        getAtom(originalRes.atoms() , "CG1"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "CD1").setXYZ(xyza[0],  xyza[1], xyza[2]);
}

 protected static void buildLYS(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
     /* Atom LYS CG CH2E */
     double bond=1.520;     /*       CH2E-CH2E */
    double angle=(Math.PI/180)*114.1; /*  CH1E-CH2E-CH2E */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG").setXYZ(xyza[0],  xyza[1], xyza[2]);

     /* Atom LYS CD CH2E */
       bond=1.520;     /*       CH2E-CH2E */
       angle=(Math.PI/180)*111.3; /*  CH2E-CH2E-CH2E */

        getAtom_xyza(xyza,  bond, angle, chi[1],
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "CD").setXYZ(xyza[0],  xyza[1], xyza[2]);

     /* Atom LYS CE CH2E */
     bond=1.520;     /*       CH2E-CH2E */
    angle=(Math.PI/180)*111.3; /*  CH2E-CH2E-CH2E */

        getAtom_xyza(xyza,  bond, angle, chi[2],
        getAtom(originalRes.atoms() , "CD"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
        
        getAtom(originalRes.atoms() , "CE").setXYZ(xyza[0],  xyza[1], xyza[2]);

     /* Atom LYS NZ NH3 */
      bond=1.489;     /*       CH2E-NH3 */
    angle=(Math.PI/180)*111.9; /*  CH2E-CH2E-NH3 */

        getAtom_xyza(xyza,  bond, angle, chi[3],
        getAtom(originalRes.atoms() , "CE"), getAtom(originalRes.atoms() , "CD"), getAtom(originalRes.atoms() , "CG"));
        
        getAtom(originalRes.atoms() , "NZ").setXYZ(xyza[0],  xyza[1], xyza[2]);
    }

    protected static void buildLEU(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
        /* Atom LEU CG CH1E */
        double bond = 1.530; /*       CH1E-CH2E */
        double angle = (Math.PI/180) * 116.3; /*  CH1E-CH2E-CH1E */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom LEU CD1 CH3E */
        bond = 1.521; /*       CH1E-CH3E */
        angle = (Math.PI/180) * 110.7; /*  CH2E-CH1E-CH3E */

        getAtom_xyza(xyza,  bond, angle, chi[1],
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "CD1").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom LEU CD2 CH3E */
        bond = 1.521; /*       CH1E-CH3E */
        angle = (Math.PI/180) * 110.7; /*  CH2E-CH1E-CH3E */

        getAtom_xyza(xyza,  bond, angle, chi[1]+(Math.PI * 2.0 / 3.0),
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "CD2").setXYZ(xyza[0],  xyza[1], xyza[2]);
    }

     protected static void buildMET(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
         /* Atom MET CG CH2E */
         double  bond=1.520;     /*       CH2E-CH2E */
         double angle=(Math.PI/180)*114.1; /*  CH1E-CH2E-CH2E */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG").setXYZ(xyza[0],  xyza[1], xyza[2]);

      /* Atom MET SD SM */
      bond=1.803;     /*       CH2E-SM */
      angle=(Math.PI/180)*112.7; /*  CH2E-CH2E-SM */

        getAtom_xyza(xyza,  bond, angle, chi[1],
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "SD").setXYZ(xyza[0],  xyza[1], xyza[2]);

      /* Atom MET CE CH3E */
      bond=1.791;     /*       SM-CH3E */
      angle=(Math.PI/180)*100.9; /*  CH2E-SM-CH3E */

        getAtom_xyza(xyza,  bond, angle, chi[2],
        getAtom(originalRes.atoms() , "SD"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
        
        getAtom(originalRes.atoms() , "CE").setXYZ(xyza[0],  xyza[1], xyza[2]);
 }

 protected static void buildASN(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
     /* Atom ASN CG C */
     double bond=1.516;     /*       CH2E-C */
     double angle=(Math.PI/180)*112.6; /*  CH1E-CH2E-C */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG").setXYZ(xyza[0],  xyza[1], xyza[2]);


     /* Atom ASN OD1 O */
     bond=1.231;     /*       C-O */
     angle=(Math.PI/180)*120.8; /*  CH2E-C-O */

       getAtom_xyza(xyza,  bond, angle, chi[1],
       getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
       
       getAtom(originalRes.atoms() , "OD1").setXYZ(xyza[0],  xyza[1], xyza[2]);

     /* Atom ASN ND2 NH2 */
     bond=1.328;     /*       C-NH2 */
     angle=(Math.PI/180)*116.4; /*  CH2E-C-NH2 */

       getAtom_xyza(xyza,  bond, angle, chi[1]-Math.PI,
       getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
       
       getAtom(originalRes.atoms() , "ND2").setXYZ(xyza[0],  xyza[1], xyza[2]);

     /* Atom ASN HD21 NH2 */
     bond=1.0;     /*       C-NH2 */
     angle=(Math.PI/180)*119; /*  CH2E-C-NH2 */

       getAtom_xyza(xyza,  bond, angle, 0,
       getAtom(originalRes.atoms() , "ND2"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
       
       getAtom(originalRes.atoms() , "HD21").setXYZ(xyza[0],  xyza[1], xyza[2]);

     /* Atom ASN HD22 NH2 */
     bond=1.0;     /*       C-NH2 */
     angle=(Math.PI/180)*119; /*  CH2E-C-NH2 */

       getAtom_xyza(xyza,  bond, angle, Math.PI,
       getAtom(originalRes.atoms() , "ND2"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
       
       getAtom(originalRes.atoms() , "HD22").setXYZ(xyza[0],  xyza[1], xyza[2]);


 }

 protected static void buildPRO(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
     /* Atom PRO CG CH2P */
     double bond=1.492;     /*       CH2E-CH2P */
     double angle=(Math.PI/180)*104.5; /*  CH1E-CH2E-CH2P */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG").setXYZ(xyza[0],  xyza[1], xyza[2]);

     /* Atom PRO CD CH2P */
     bond=1.503;     /*       CH2P-CH2P */
    angle=(Math.PI/180)*106.1; /*  CH1E-CH2P-CH2P */

        getAtom_xyza(xyza,  bond, angle, chi[1],
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "CD").setXYZ(xyza[0],  xyza[1], xyza[2]);
  }

  protected static void buildGLN(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
      /* Atom GLN CG CH2E */
      double bond=1.520;     /*       CH2E-CH2E */
      double angle=(Math.PI/180)*114.1; /*  CH1E-CH2E-CH2E */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG").setXYZ(xyza[0],  xyza[1], xyza[2]);


      /* Atom GLU CD C */
      bond= 1.516;     /*       CH2E-C */
      angle=(Math.PI/180)*112.6; /*  CH2E-CH2E-C */

        getAtom_xyza(xyza,  bond, angle, chi[1],
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "CD").setXYZ(xyza[0],  xyza[1], xyza[2]);

      /* Atom GLN OE1 O */
      bond=1.231;     /*       C-O */
    angle=(Math.PI/180)*120.8; /*  CH2E-C-O */

       getAtom_xyza(xyza,  bond, angle, chi[2],
       getAtom(originalRes.atoms() , "CD"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
       
       getAtom(originalRes.atoms() , "OE1").setXYZ(xyza[0],  xyza[1], xyza[2]);

      /* Atom GLN NE2 NH2 */
      bond=1.328;     /*       C-NH2 */
      angle=(Math.PI/180)*116.4; /*  CH2E-C-NH2 */

       getAtom_xyza(xyza,  bond, angle, chi[2]-Math.PI,
       getAtom(originalRes.atoms() , "CD"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
       
       getAtom(originalRes.atoms() , "NE2").setXYZ(xyza[0],  xyza[1], xyza[2]);

      /* Atom GLN NE21 NH2 */
      bond=1.0;     /*       C-NH2 */
      angle=(Math.PI/180)*119.0; /*  CH2E-C-NH2 */

       getAtom_xyza(xyza,  bond, angle, 0,
       getAtom(originalRes.atoms() , "NE2"), getAtom(originalRes.atoms() , "CD"), getAtom(originalRes.atoms() , "CG"));
       
       getAtom(originalRes.atoms() , "HE21").setXYZ(xyza[0],  xyza[1], xyza[2]);

      /* Atom GLN NE22 NH2 */
      bond=1.0;     /*       C-NH2 */
      angle=(Math.PI/180)*119.0; /*  CH2E-C-NH2 */

       getAtom_xyza(xyza,  bond, angle, Math.PI,
       getAtom(originalRes.atoms() , "NE2"), getAtom(originalRes.atoms() , "CD"), getAtom(originalRes.atoms() , "CG"));
       
       getAtom(originalRes.atoms() , "HE22").setXYZ(xyza[0],  xyza[1], xyza[2]);
    }



    protected static void buildARG(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
        /* Atom LYS CG CH2E */
        double bond=1.520;     /*       CH2E-CH2E */
       double angle=(Math.PI/180)*114.1; /*  CH1E-CH2E-CH2E */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG").setXYZ(xyza[0],  xyza[1], xyza[2]);


        /* Atom LYS CD CH2E */
          bond=1.520;     /*       CH2E-CH2E */
          angle=(Math.PI/180)*111.3; /*  CH2E-CH2E-CH2E */

        getAtom_xyza(xyza,  bond, angle, chi[1],
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "CD").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom ARG NE NH1 */
        bond=1.461;     /*      CH2E-NH1 */
        angle=(Math.PI/180)*112.0; /* CH2E-CH2E-NH1 */

        getAtom_xyza(xyza,  bond, angle, chi[2],
        getAtom(originalRes.atoms() , "CD"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
        
        getAtom(originalRes.atoms() , "NE").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom ARG HE */
        bond=1.0;     /*      CH2E-NH1 */
        angle=(Math.PI/180)*117.58; /* CH2E-CH2E-NH1 */

        getAtom_xyza(xyza,  bond, angle, chi[3]-Math.PI,
        getAtom(originalRes.atoms() , "NE"), getAtom(originalRes.atoms() , "CD"), getAtom(originalRes.atoms() , "CG"));
        
        getAtom(originalRes.atoms() , "HE").setXYZ(xyza[0],  xyza[1], xyza[2]);


        /* Atom ARG CZ C */
        bond=1.329;     /*       NH1-C */
        angle=(Math.PI/180)*124.2; /*  CH2E-NH1-C */

        getAtom_xyza(xyza,  bond, angle, chi[3],
        getAtom(originalRes.atoms() , "NE"), getAtom(originalRes.atoms() , "CD"), getAtom(originalRes.atoms() , "CG"));
        
        getAtom(originalRes.atoms() , "CZ").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom ARG NH1 NC2 */
        bond=1.326;     /*      C-NC2 */
        angle=(Math.PI/180)*120.0; /*  NH1-C-NC2 */

        getAtom_xyza(xyza,  bond, angle, 0,
        getAtom(originalRes.atoms() , "CZ"), getAtom(originalRes.atoms() , "NE"), getAtom(originalRes.atoms() , "CD"));
        
        getAtom(originalRes.atoms() , "NH1").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom ARG NH2 NH2 */
        bond=1.326;     /*      C-NC2 */
        angle=(Math.PI/180)*120.0; /*  NH1-C-NC2 */

        getAtom_xyza(xyza,  bond, angle, Math.PI,
        getAtom(originalRes.atoms() , "CZ"), getAtom(originalRes.atoms() , "NE"), getAtom(originalRes.atoms() , "CD"));
        
        getAtom(originalRes.atoms() , "NH2").setXYZ(xyza[0],  xyza[1], xyza[2]);
    }

    protected static void buildSER(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
        /* Atom SER OG OH1 */
        double bond=1.417;     /*       CH2E-OH1 */
        double angle=(Math.PI/180)*111.1; /*  CH1E-CH2E-OH1 */


        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "OG").setXYZ(xyza[0],  xyza[1], xyza[2]);
    }

    protected static void buildTHR(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
        /* Atom THR OG1 OH1 */
        double bond=1.433;     /*       CH1E-OH1 */
        double angle=(Math.PI/180)*109.6; /*  CH1E-CH1E-OH1 */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "OG1").setXYZ(xyza[0],  xyza[1], xyza[2]);
        /* Atom THR CG2 CH3E */
        bond=1.521;     /*       CH1E-CH3E */
        angle=(Math.PI/180)*110.5; /*  CH1E-CH1E-CH3E */

        getAtom_xyza(xyza,  bond, angle, chi[0]-(Math.PI*2/3),
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG2").setXYZ(xyza[0],  xyza[1], xyza[2]);
    }

    protected static void buildVAL(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
        /* Atom VAL CG1 CH3E */
        double bond=1.521;     /*       CH1E-CH3E */
        double angle=(Math.PI/180)*110.5; /*  CH1E-CH1E-CH3E */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG1").setXYZ(xyza[0],  xyza[1], xyza[2]);
        /* Atom VAL CG2 CH3E */
        bond=1.521;     /*       CH1E-CH3E */
        angle=(Math.PI/180)*110.5; /*  CH1E-CH1E-CH3E */

        getAtom_xyza(xyza,  bond, angle, chi[0]+(Math.PI*2/3),
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG2").setXYZ(xyza[0],  xyza[1], xyza[2]);
    }

    protected static void buildTRP(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
        /* Atom TRP CG C5W */
        double bond=1.498;     /*       CH2E-C5W */
        double angle=(Math.PI/180)*113.6; /*  CH1E-CH2E-C5W */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom TRP CD1 CR1E */
         bond=1.365;     /*       C5W-CR1E */
         angle=(Math.PI/180)*126.934; /*  CH2E-C5W-CR1E */

        getAtom_xyza(xyza,  bond, angle, chi[1],
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "CD1").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom TRP CD2 CW */
         bond=1.433;     /*       C5W-CW */
         angle=(Math.PI/180)*126.578; /*  CH2E-C5W-CW */

        getAtom_xyza(xyza,  bond, angle, chi[1]-Math.PI,
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "CD2").setXYZ(xyza[0],  xyza[1], xyza[2]);

         /* Atom TRP NE1 NH1 */
         bond=1.374;     /*      CR1E-NH1 */
         angle=(Math.PI/180)*110.2; /*  C5W-CR1E-NH1 */

        getAtom_xyza(xyza,  bond, angle, Math.PI,
        getAtom(originalRes.atoms() , "CD1"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
        
        getAtom(originalRes.atoms() , "NE1").setXYZ(xyza[0],  xyza[1], xyza[2]);

         /* Atom TRP NE1 NH1 */
         bond=1.0;     /*      CR1E-NH1 */
         angle=(Math.PI/180)*125.2; /*  C5W-CR1E-NH1 */

        getAtom_xyza(xyza,  bond, angle, Math.PI,
        getAtom(originalRes.atoms() , "NE1"), getAtom(originalRes.atoms() , "CD1"), getAtom(originalRes.atoms() , "CG"));
        
        getAtom(originalRes.atoms() , "HE1").setXYZ(xyza[0],  xyza[1], xyza[2]);


         /* Atom TRP CE2 CW */
         bond=1.4114;     /*      CW-CW */
         angle=(Math.PI/180)*107.1; /*  C5W-CW-CW */

        getAtom_xyza(xyza,  bond, angle, Math.PI,
        getAtom(originalRes.atoms() , "CD2"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
        
        getAtom(originalRes.atoms() , "CE2").setXYZ(xyza[0],  xyza[1], xyza[2]);

         /* Atom TRP CE3 CR1E */
         bond=1.401;     /*     CW-CR1E */
         angle=(Math.PI/180)*133.885; /*  CW-CW-CR1E */

        getAtom_xyza(xyza,  bond, angle, 0,
        getAtom(originalRes.atoms() , "CD2"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
        
        getAtom(originalRes.atoms() , "CE3").setXYZ(xyza[0],  xyza[1], xyza[2]);

         /* Atom TRP CZ2 CR1W */
         bond=1.3982;     /*     CW-CR1W */
         angle=(Math.PI/180)*122.332; /*  CW-CW-CR1W */

        getAtom_xyza(xyza,  bond, angle, Math.PI,
        getAtom(originalRes.atoms() , "CE2"), getAtom(originalRes.atoms() , "CD2"), getAtom(originalRes.atoms() , "CG"));
        
        getAtom(originalRes.atoms() , "CZ2").setXYZ(xyza[0],  xyza[1], xyza[2]);

         /* Atom TRP CZ3 CR1E */
         bond=1.391;     /*     CR1E-CR1E */
         angle=(Math.PI/180)*118.654; /*  CW-CR1E-CR1E */

        getAtom_xyza(xyza,  bond, angle, Math.PI,
        getAtom(originalRes.atoms() , "CE3"), getAtom(originalRes.atoms() , "CD2"), getAtom(originalRes.atoms() , "CG"));
        
        getAtom(originalRes.atoms() , "CZ3").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom TRP CH2 CR1W */
        bond=1.368;     /*     CR1W-CR1W */
        angle=(Math.PI/180)*117.5; /*  CW-CR1W-CR1W */

        getAtom_xyza(xyza,  bond, angle, 0,
        getAtom(originalRes.atoms() , "CZ3"), getAtom(originalRes.atoms() , "CE3"), getAtom(originalRes.atoms() , "CD2"));
        
        getAtom(originalRes.atoms() , "CH2").setXYZ(xyza[0],  xyza[1], xyza[2]);
    }

    protected static void buildTYR(Residue originalRes,double[] chi) {
    	buildALA(originalRes);
        /* Atom TYR CG CY */
        double bond=1.512;     /*       CH2E-CY */
        double angle=(Math.PI/180)*113.9; /*  CH1E-CH2E-CY */

        getAtom_xyza(xyza,  bond, angle, chi[0],
        getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CG").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom TYR CD1 CR1E */
        bond=1.389;     /*       CY-CR1E */
        angle=(Math.PI/180)*120.8; /*  CH2E-CY-CR1E */

        getAtom_xyza(xyza,  bond, angle, chi[1],
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "CD1").setXYZ(xyza[0],  xyza[1], xyza[2]);

        /* Atom TYR CD2 CR1E */
        bond=1.389;     /*      CY-CR1E */
        angle=(Math.PI/180)*120.8; /* CH2E-CY-CR1E */

        getAtom_xyza(xyza,  bond, angle, chi[1]-Math.PI,
        getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"), getAtom(originalRes.atoms() , "CA"));
        
        getAtom(originalRes.atoms() , "CD2").setXYZ(xyza[0],  xyza[1], xyza[2]);

         /* Atom TYR CE1 CR1E */
         bond=1.382;     /*     CR1E-CR1E */
         angle=(Math.PI/180)*121.2; /*  CY-CR1E-CR1E */

        getAtom_xyza(xyza,  bond, angle, Math.PI,
        getAtom(originalRes.atoms() , "CD1"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
        
        getAtom(originalRes.atoms() , "CE1").setXYZ(xyza[0],  xyza[1], xyza[2]);

         /* Atom TYR CE2 CR1E */
         bond=1.382;     /*     CR1E-CR1E */
         angle=(Math.PI/180)*121.2; /*  CY-CR1E-CR1E */

        getAtom_xyza(xyza,  bond, angle, Math.PI,
        getAtom(originalRes.atoms() , "CD2"), getAtom(originalRes.atoms() , "CG"), getAtom(originalRes.atoms() , "CB"));
        
        getAtom(originalRes.atoms() , "CE2").setXYZ(xyza[0],  xyza[1], xyza[2]);

         /* Atom TYR CZ CY2 */
         bond=1.378;     /*       CR1E-CY2 */
         angle=(Math.PI/180)*119.6; /*  CR1E-CR1E-CY2 */

        getAtom_xyza(xyza,  bond, angle, 0,
        getAtom(originalRes.atoms() , "CE2"), getAtom(originalRes.atoms() , "CD2"), getAtom(originalRes.atoms() , "CG"));
        
        getAtom(originalRes.atoms() , "CZ").setXYZ(xyza[0],  xyza[1], xyza[2]);

         /* Atom TYR OH OH1 */
         bond=1.376;     /*       CY2-OH1 */
         angle=(Math.PI/180)*119.9; /*  CR1E-CY2-OH1 */

        getAtom_xyza(xyza,  bond, angle, Math.PI,
        getAtom(originalRes.atoms() , "CZ"), getAtom(originalRes.atoms() , "CE1"), getAtom(originalRes.atoms() , "CD1"));
        
        getAtom(originalRes.atoms() , "OH").setXYZ(xyza[0],  xyza[1], xyza[2]);
     }

/**
 * This will put the CB of the residue in the coordinates of the side-chain centroid.
 * @param originalRes
 */
    public static void buildCentroid(Residue originalRes) {
    	// Build CB
        double bond=0;     
        double angle=0;
        double tor=0;
    

        switch (originalRes.type) {
        case 0:  bond = 1.522661933747001; //A
        angle = 1.9209171643112692;
        tor = 2.1338779181999596;
        break;
        case 1:  bond = 1.9301945789274642; //C
        angle = 2.0516803398041326;
        tor = 2.228430522427323;
        break;
        case 2:  bond = 2.175733745355916; //D
        angle = 2.116376516980476;
        tor = 2.3110072204242287;
        break;
        case 3:  bond = 2.745200734723785; //E
        angle = 2.1280907062775167;
        tor = 2.3599671591898126;
        break;
        case 4:  bond = 2.726200320738368; //F
        angle = 2.1488101592255076;
        tor = 2.448274892623533;
        break;
        case 5:  return ;//G
        case 6:  bond = 2.567985094096986; //H
        angle = 2.1943006994808583;
        tor = 2.3500027269561996;
        break;
        case 7:  bond = 2.27664187151238; //I
        angle = 2.0966566204119315;
        tor = 2.2803008321709832;
        break;
        case 8:  bond = 3.199690207299597; //K
        angle = 2.1464050693298535;
        tor = 2.3612772698633764;
        break;
        case 9:  bond = 2.494092206852833; //L
        angle = 2.1740945329611927;
        tor = 2.4463879018234946;
        break;
        case 10:  bond = 2.711514408220479; //M
        angle = 2.1195853435121395;
        tor = 2.431271717014755;
        break;
        case 11:  bond = 2.1665234701808926; //N
        angle = 2.1479162129543363;
        tor = 2.3310699215193735;
        break;
        case 12:  bond = 1.8416736742860633; //P
        angle = 2.1799102731169517;
        tor = 1.2589175769756424;
        break;
        case 13:  bond = 2.7311910723192954; //Q
        angle = 2.140339921228663;
        tor = 2.4031147069070062;
        break;
        case 14:  bond = 3.5842923428997655; //R
        angle = 2.100086132765975;
        tor = 2.337839390338103;
        break;
        case 15:  bond = 1.7845213050435664; //S
        angle = 1.8638353477337217;
        tor = 2.082561926825327;
        break;
        case 16:  bond = 1.8945986458317374; //T
        angle = 2.002610648669302;
        tor = 2.1210569564568473;
        break;
        case 17:  bond = 1.933950772000282; //V
        angle = 1.9946619685057034;
        tor = 2.291833818501636;
        break;
        case 18:  bond = 2.9804216764794638; //W
        angle = 2.0431735676378757;
        tor = 2.5252664877304194;
        break;
        case 19:  bond = 2.8931077870708184; //Y
        angle = 2.1404097084546954;
        tor = 2.437885087976493;
        break;
        default: throw new RuntimeException("wrong residue type: "+originalRes.number + " of type " + originalRes.type);
        }
    	
    	
        getAtom_xyza(xyza,  bond, angle, tor,
        getAtom(originalRes.atoms() , "CA"), getAtom(originalRes.atoms() , "C"), getAtom(originalRes.atoms() , "N"));
        
        getAtom(originalRes.atoms() , "CB").setXYZ(xyza[0],  xyza[1], xyza[2]);    	
    }

    
    /**
     * Returns the specified atom from a residue.
     **/
    private static Atom getAtom(AtomList al, String atomName){
        Iterator atoms;
        Atom atom;
        atoms = al.iterator();
        while ((atom = (Atom) atoms.next()) != null)
           if (atom.name.equals(atomName)) {
              return atom;
           }
//        System.out.println("\n\n\n\n" + al);
//        System.out.println("\n\nWarning!! Could not find atom name: " + atomName + " in the list above");
        return null;
    }
    
    

}

