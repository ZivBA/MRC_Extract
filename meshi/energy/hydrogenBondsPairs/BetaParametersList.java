package meshi.energy.hydrogenBondsPairs;

import meshi.util.filters.Filter;

/**
 * Created by IntelliJ IDEA.
 * User: amilev
 * Date: 06/02/2006
 * Time: 14:04:55
 * To change this template use File | Settings | File Templates.
 */
public class BetaParametersList extends HydrogenBondsPairsParametersList {



    //------------------------------------ constructors ----------------------------
      public BetaParametersList(){
          super(new IsHydrogenBondsPairsParameters());
     }

     public BetaParametersList(String parametersFileName) {
        super(parametersFileName);
    }

    //---------------------------------------------------------------------------------------------
    private static class IsHydrogenBondsPairsParameters implements Filter {
    public boolean accept(Object obj) {
        return (obj instanceof HydrogenBondsPairsParameters);
    }
    }
}
