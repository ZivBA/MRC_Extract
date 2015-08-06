package meshi.util.crossLinking;

import meshi.parameters.AtomTypes;

public interface ResidueMasses extends AtomTypes {
    public static final double[] MW_AA = 
       {/*Ala*/	71.037114,
    	/*Cys*/	174.04628, // 103.00918 for un-alkylated cys , 161.014680 for the one methyl less
    	/*Asp*/	115.02694,
    	/*Glu*/	129.04259,
    	/*Phe*/	147.06841,
    	/*Gly*/	57.021464,
    	/*His*/	137.05891,
    	/*Ile*/	113.08406,
    	/*Lys*/	128.09496,
    	/*Leu*/	113.08406,
    	/*Met*/	131.04048,
    	/*Asn*/	114.04293,
    	/*Pro*/	97.052764,
    	/*Gln*/	128.05858,
    	/*Arg*/	156.10111,
    	/*Ser*/	87.032028,
    	/*Thr*/	101.04768,
    	/*Val*/	99.068414,
    	/*Trp*/	186.07931,
    	/*Tyr*/	163.06333};
    
    public static final double MW_H = 1.0078250;
    public static final double MW_O = 15.994915;
    public static final double MW_OH = 17.00274;
    public static final double MW_H2O = 18.01057;
    public static final double MW_CL = 138.0681; // BS3
    
    public static final int atomLinked1 = KNZ;
    public static final int atomLinked2 = KNZ;

    public static final double[][] MW_AA_MODIfICATIONS = 
    {/*Ala*/ {71.037114},
 	/*Cys*/ {174.0463}, // 103.00918 for un-alkylated cys , 160.030648 ,161.051049, 174.0463, 175.0667
 	/*Asp*/	{115.02694},
 	/*Glu*/	{129.04259},
 	/*Phe*/	{147.06841},
 	/*Gly*/	{57.021464},
 	/*His*/	{137.05891},
 	/*Ile*/	{113.08406},
 	/*Lys*/	{128.09496/*, 284.1736042, 387.23691*/},// The mono-link modifications are: 284.1736042 (reaction with water) and 387.23691 (reaction with Tris)
 	/*Leu*/	{113.08406},
 	/*Met*/	{131.04048}, // Met oxidation is: 147.03540
 	/*Asn*/	{114.04293},
 	/*Pro*/	{97.052764},
 	/*Gln*/	{128.05858},
 	/*Arg*/	{156.10111},
 	/*Ser*/	{87.032028},
 	/*Thr*/	{101.04768},
 	/*Val*/	{99.068414},
 	/*Trp*/	{186.07931},
 	/*Tyr*/	{163.06333},
 	/*N-term*/ {0/*,156.0786442,259.14195*/}};  // The mono-link modifications are: 156.0786442 (reaction with water) and 259.14195 (reaction with Tris) 

    
    
}
