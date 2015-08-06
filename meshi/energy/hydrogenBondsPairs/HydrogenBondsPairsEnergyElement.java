package meshi.energy.hydrogenBondsPairs;

import meshi.energy.NonBondedEnergyElement;
import meshi.energy.hydrogenBond.HB_DistanceAttribute;
import meshi.energy.hydrogenBond.HydrogenBondsEnergyElement;
import meshi.geometry.Distance;
import meshi.molecularElements.Atom;

public class HydrogenBondsPairsEnergyElement extends NonBondedEnergyElement {

    //----------------------------------------- data ------------------------------------

    BetaParametersList betaParametersList;
    HelixParametersList helixParametersList;
    double punish,hpunish;
    private Atom oAtom1, hAtom1, oAtom2, hAtom2;    
    private double weight;    
    Distance HOelement1,HOelement2;
    HB_DistanceAttribute distanceAttributes1,distanceAttributes2;
    private double deDxOAtom1, deDyOAtom1, deDzOAtom1,deDxOAtom2, deDyOAtom2, deDzOAtom2;
    private double deDxHAtom1, deDyHAtom1, deDzHAtom1,deDxHAtom2, deDyHAtom2, deDzHAtom2;    
    private double energy = 0;

    private boolean free = true;
    private double factor; //when factor >0 this patern has been seen in the datebase
    private int[] specialDisArray = null;
    private boolean antiParallel = true;
    //when factor <0 this patern is rare in the datebase


    public final Atom oAtom1(){return oAtom1;}
        
    //----------------------------------------- constructors --------------------------------
    
    public HydrogenBondsPairsEnergyElement() {}
    public HydrogenBondsPairsEnergyElement(HelixParametersList helixParametersList,
                                           BetaParametersList betaParametersList,
                                           double weight,
                                           double punish,
                                           double hpunish){
        this.helixParametersList= helixParametersList;
        this.betaParametersList = betaParametersList;
        this.weight = weight;
        this.punish = punish;
        this.hpunish = hpunish;
    }

    public HydrogenBondsPairsEnergyElement(HelixParametersList helixParametersList,
                                           BetaParametersList betaParametersList,
                                           double weight,
                                           double punish,
                                           double hpunish,
                                           int[] specialDisArray,
                                           boolean antiParallel){
        this(helixParametersList,betaParametersList,weight,punish,hpunish);
        this.specialDisArray = specialDisArray ;
        this.antiParallel = antiParallel;
    }


    //------------------------------------------ methods -------------------------------------
    
    public void set(Object obj){//obj should be PairOfHydrogenBondsElements
        free = false;
        PairOfHydrogenBondsElements pair = (PairOfHydrogenBondsElements) obj;
        if(! pair.lookAtThisPair)
            throw new RuntimeException("problem in HydrogenBondsPairsEnergyElement:this pair should not be looked at: "+pair);        

        atoms = pair.atoms();
        this.HOelement1 = pair.HOelement1;
        this.HOelement2 = pair.HOelement2;
        distanceAttributes1 = (HB_DistanceAttribute)(HOelement1.getAttribute(HB_DistanceAttribute.key));
        distanceAttributes2 = (HB_DistanceAttribute)(HOelement2.getAttribute(HB_DistanceAttribute.key));
        if (distanceAttributes1 == null || distanceAttributes2 == null)
                   throw new RuntimeException("problem in HydrogenBondsPairsEnergyElement: one of the element is null");
        oAtom1 = distanceAttributes1.oAtom;
        hAtom1 = distanceAttributes1.hAtom;
        oAtom2 = distanceAttributes2.oAtom;
        hAtom2 = distanceAttributes2.hAtom;
        //System.out.println(oAtom1);
        HydrogenBondsPairsParameters parameters = null;
        if(pair.isHelixPair() ) {
           // System.out.println("HydrogenBondsPairsEnergyElement: helixPair");
            parameters = (HydrogenBondsPairsParameters)helixParametersList.parameters(pair);
        }
        else if(pair.isBetaPair() )
          parameters = (HydrogenBondsPairsParameters)betaParametersList.parameters(pair ) ;//  System.out.println("HydrogenBondsEnergyElement: betaPair");
         else if(pair. isOneBetta() )
            System.out.println("HydrogenBondsEnergyElement: oneBeta");
         else if(pair.isOneHelix() )
            System.out.println("HydrogenBondsEnergyElement: oneHelix");
        else
            throw new RuntimeException("pair should be beta/helix pair");
        //System.out.println("hydrogenBondsPairsEnergyElement:updateEnergy:parameters: /n "+parameters);
        double parameterValue = parameters.value();
        pair.setPairValue(parameterValue );
        int hhDistance = parameters.h1h2SeqGap;
        int ooDistance = parameters .o1o2SeqGap;
        if (parameters.firstHBSeqGap == -4 &
                parameters .secondHBSeqGap == -4 &
                parameters .h1h2SeqGap == 1 &
                parameters .o1o2SeqGap == 1 &
                parameters .h1o2SeqGap == 5 &
                parameters .o1h2SeqGap == -3)  //TODO check it
                factor = 1000;             //This is a helix !
           else if(ooDistance ==  0){                  //TODO add a spacial punish to this case, if needed (if results indicates that there are more then expected o atoms with more the one hbonds
               factor = punish; //TODO check         
            }
            else if (hhDistance == 0){
                if(pair.sheetPair(HOelement1 ) | pair .sheetPair(HOelement2 ) )      //TODO what if it's not a pair sheet ?
                      factor = hpunish;//*10;
             else factor = hpunish;          //negative
        }
        else if (parameterValue < 0){   //TODO should it be < 1  or <= 1 as it was
            if (ooDistance ==  0)        
                      System.out .println("ooDistance ==  0 :"+pair);
            factor = punish; //negative
        }
         else if(specialDisArray != null){
            if(antiParallel){
                if (!antiParallel(parameters ))
                    factor = 0;
                else if(antiParallelshift(hAtom1 .residueNumber(),oAtom1 .residueNumber()))
                    factor = -1;
                else
                    factor = parameterValue/10;
            }
            else
                if( ! parallel(parameters ) )
                    factor = 0;
                else if(parallelshift(hAtom1 .residueNumber(),oAtom1 .residueNumber()))
                    factor = -1;
                else
                factor = parameterValue/10;
        }
        else
            factor = parameterValue/10;//positive          //TODO check different scales.
        if(distanceAttributes1 != HOelement1.getAttribute(HB_DistanceAttribute.key))
            System.out.println("this is wird");
        pair.setPairValue(factor);
      }

    private boolean antiParallel(HydrogenBondsPairsParameters parameters){
        return (parameters .h1h2SeqGap ==-1*parameters.o1o2SeqGap & parameters.h1o2SeqGap == -1*parameters.o1h2SeqGap) ;
    }

    private boolean parallel(HydrogenBondsPairsParameters parameters){
        return  ((parameters .h1h2SeqGap == parameters.o1o2SeqGap & (parameters .h1h2SeqGap == 2 | parameters .h1h2SeqGap == 4) & parameters.h1o2SeqGap == -1*parameters.o1h2SeqGap) |
                          (parameters .h1h2SeqGap == -1*parameters.o1o2SeqGap & parameters.h1o2SeqGap == parameters.o1h2SeqGap +2  & (parameters.h1o2SeqGap == 4 | parameters.h1o2SeqGap == 2 | parameters.h1o2SeqGap == 0 | parameters.h1o2SeqGap == -2)));

    }

    private boolean antiParallelshift (int resNum1,int resNum2){
        return Math.abs(resNum1 - resNum2) % 2 != Math.abs(specialDisArray[0] - specialDisArray[1]) % 2;
    }

     private boolean parallelshift (int resNum1,int resNum2){
        return Math.abs(resNum1 - resNum2) % 2 != Math.abs(specialDisArray[0] - (specialDisArray[1]+1)) % 2;
    }

    public void freeElenet(){
        if(distanceAttributes1 != HOelement1.getAttribute(HB_DistanceAttribute.key))
            System.out.println("this is wird");
        free = true;
        atoms = null;
        HOelement1 = HOelement2 = null;
        distanceAttributes2 = distanceAttributes1 = null;
        factor = 1/0.0;
        deDxHAtom1 = deDxHAtom2 =deDxOAtom1 =deDxOAtom2 = 1/0.0;
        deDyHAtom1 = deDyHAtom2 =deDyOAtom1 =deDyOAtom2 =1/0.0;
        deDzHAtom1 =deDzHAtom2 =deDzOAtom1 =deDzOAtom2 =1/0.0;
        energy =1/0.0;
        oAtom1 =oAtom2 =hAtom1 =hAtom2=null;
    }

    /* (non-Javadoc)
    * @see meshi.energy.EnergyElement#evaluate()
    */
    public double evaluate() {
       // int   x;
      //  if(hAtom1.number() ==29 | hAtom2.number() == 29)
      //      x = 1;
        if(distanceAttributes1 != HOelement1.getAttribute(HB_DistanceAttribute.key))
            System.out.println("this is wird");
        updateEnergy();
        updateAtoms();
        if(distanceAttributes1 != HOelement1.getAttribute(HB_DistanceAttribute.key))
            System.out.println("this is wird");
        return energy * weight;
    }
	
    public double updateEnergy() {

        HydrogenBondsEnergyElement hbElement = new HydrogenBondsEnergyElement();
        hbElement.set(HOelement1);
        if(distanceAttributes1 != HOelement1.getAttribute(HB_DistanceAttribute.key))
            System.out.println("this is wird");
        if(distanceAttributes1 != hbElement.hb_Attribute )
            System.out.println("this is wird");
        hbElement.updateEnergy();
         hbElement.freeElement();
        hbElement.set(HOelement2); //TODO try use new element
         hbElement.updateEnergy();
         hbElement.freeElement();
        //TODO take care of this in more elegant way, now there are double calculations of HB energy

        double e1, e2;
        double dE1DxOAtom1, dE1DyOAtom1, dE1DzOAtom1;//, dE1DxOAtom2, dE1DyOAtom2, dE1DzOAtom2;
        double dE2DxOAtom2, dE2DyOAtom2, dE2DzOAtom2;//, dE2DxOAtom1, dE2DyOAtom1, dE2DzOAtom1;
        double dE1DxHAtom1, dE1DyHAtom1, dE1DzHAtom1;//, dE1DxHAtom2, dE1DyHAtom2, dE1DzHAtom2;
        double dE2DxHAtom2, dE2DyHAtom2, dE2DzHAtom2;//, dE2DxHAtom1, dE2DyHAtom1, dE2DzHAtom1;
            
        e1 = distanceAttributes1.energy;        
        
        dE1DxOAtom1 = distanceAttributes1.deDxOAtom;
        dE1DyOAtom1 = distanceAttributes1.deDyOAtom;        
        dE1DzOAtom1 = distanceAttributes1.deDzOAtom;
        dE1DxHAtom1 = distanceAttributes1.deDxHAtom;
        dE1DyHAtom1 = distanceAttributes1.deDyHAtom;
        dE1DzHAtom1 = distanceAttributes1.deDzHAtom;
  
        e2 = distanceAttributes2.energy;                

        dE2DxOAtom2 = distanceAttributes2.deDxOAtom;
        dE2DyOAtom2 = distanceAttributes2.deDyOAtom;
        dE2DzOAtom2 = distanceAttributes2.deDzOAtom;
        dE2DxHAtom2 = distanceAttributes2.deDxHAtom;
        dE2DyHAtom2 = distanceAttributes2.deDyHAtom;
        dE2DzHAtom2 = distanceAttributes2.deDzHAtom;
        deDxOAtom1 =  -1 * factor * (e2 * dE1DxOAtom1);
        deDyOAtom1 =  -1 * factor * (e2 * dE1DyOAtom1);
        deDzOAtom1 =  -1 * factor * (e2 * dE1DzOAtom1);
        deDxHAtom1 =  -1 * factor * (e2 * dE1DxHAtom1);
        deDyHAtom1 =  -1 * factor * (e2 * dE1DyHAtom1);
        deDzHAtom1 =  -1 * factor * (e2 * dE1DzHAtom1);
                
        deDxOAtom2 =  -1 * factor * (e1 * dE2DxOAtom2 );
        deDyOAtom2 =  -1 * factor * (e1 * dE2DyOAtom2 );
        deDzOAtom2 =  -1 * factor * (e1 * dE2DzOAtom2 );
        deDxHAtom2 =  -1 * factor * (e1 * dE2DxHAtom2 );
        deDyHAtom2 =  -1 * factor * (e1 * dE2DyHAtom2 );
        deDzHAtom2 =  -1 * factor * (e1 * dE2DzHAtom2 );

        // e1, e2 are negative values => positive factor meens negative (good) energy
        //                               negative factor meens positive (bad) energy
        //System.out.println("hydrogenBondsPairsEnergyElement:updateEnergy:e1*e2 : "+e1*e2+" factor: "+factor);
        energy = -1 * factor * e1 * e2;
        if (factor < 0 && energy <0)
            System.out.println("energy: "+energy+" e1:  "+e1+" e2: "+e2+" weight: "+weight+" factor: "+factor);
        if (energy < -10000 || energy > 1000)
            System.out.println("energy: "+energy+" e1:  "+e1+" e2: "+e2+" weight: "+weight+" factor: "+factor);
        return energy;  
    }	
	
    public void updateAtoms(){
        if (! oAtom1.frozen()) {
            oAtom1.addToFx(-1 * deDxOAtom1 * weight); // force = -derivative
            oAtom1.addToFy(-1 * deDyOAtom1 * weight); // force = -derivative
            oAtom1.addToFz(-1 * deDzOAtom1 * weight); // force = -derivative       
        }
        if (! hAtom1.frozen()) {
            hAtom1.addToFx(-1 * deDxHAtom1 * weight);
            hAtom1.addToFy(-1 * deDyHAtom1 * weight);
            hAtom1.addToFz(-1 * deDzHAtom1 * weight);
        }
        if (! oAtom2.frozen()) {
            oAtom2.addToFx(-1 * deDxOAtom2 * weight); // force = -derivative
            oAtom2.addToFy(-1 * deDyOAtom2 * weight); // force = -derivative
            oAtom2.addToFz(-1 * deDzOAtom2 * weight); // force = -derivative            
        }
        if (! hAtom2.frozen()) {
            hAtom2.addToFx(-1 * deDxHAtom2 * weight);
            hAtom2.addToFy(-1 * deDyHAtom2 * weight);
            hAtom2.addToFz(-1 * deDzHAtom2 * weight);
        }
    }
	
    /* (non-Javadoc)
     * @see meshi.energy.EnergyElement#setAtoms()
     */
    protected void setAtoms() {
        throw new UnsupportedOperationException("setAtoms() may not be used by HydrogenBondsEnergyElement for "+
                                   "efficiency.") ;
    }
	
    public String toString() {
        if(free) return "HydrogenBondsPairsEnergyElement: is free (no values)";
        return "HydrogenBondsPairsEnergyElement: oAtom1 "+oAtom1.residueNumber()+" HAtom1 " +hAtom1.residueNumber()+" oAtom2 "+oAtom2.residueNumber()+" HAtom2 " +hAtom2.residueNumber();
        //return "HydrogenBondsPairsEnergyElement: factor "+factor+" "+pair;
    }		
}
