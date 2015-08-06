package meshi.energy.excludedVol;
import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.filters.Filter;

public class ExcludedVolCreator extends EnergyCreator  implements KeyWords {

    private double Rfac;
    private Filter filter = null;
 
    public ExcludedVolCreator(double weight,double Rfac) {
  	super(weight);
  	this.Rfac = Rfac;
    }
    public ExcludedVolCreator(double Rfac) {
  	super(EXCLUDED_VOL);
  	this.Rfac = Rfac;
    }

    public ExcludedVolCreator() {
  	super(EXCLUDED_VOL);
  	this.Rfac = 1.0;
    }
       
    public ExcludedVolCreator(Filter filter,double Rfac) {
	super(EXCLUDED_VOL);
	this.Rfac = Rfac;
	this.filter = filter;
    }
       

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix, 
					  CommandList commands) {
	if (parametersList== null)
	    parametersList = new ExcludedVolParametersList(parametersDirectory(commands)+
							    "/"+EXCLUDED_VOL_PARAMETERS);
	return new ExcludedVol(distanceMatrix, ( ExcludedVolParametersList) parametersList, Rfac, weight(),filter);
    }
}
