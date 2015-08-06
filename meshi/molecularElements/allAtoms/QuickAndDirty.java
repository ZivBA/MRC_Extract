package meshi.molecularElements.allAtoms;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.ResidueBuilder;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionList;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.CommandList;
import meshi.util.KeyWords;


public class QuickAndDirty implements Residues, AtomTypes , KeyWords{ /**
									 * The implemented
									 * interfaces defines the 
									 * names of atom and residue 
									 * types. 
									 **/
    private DunbrackLib lib;
    private Protein protein; 
    private DistanceMatrix distanceMatrix;
    

    public QuickAndDirty(CommandList commands) {
	lib = new DunbrackLib(commands, 0.99,2);	
	protein = null;
	distanceMatrix = null;
    }
    
    public void model(Protein protein)  {
	model(protein, protein.residues());
    }

    public void model(Protein protein, ResidueList residuesToModel)  {
    this.protein = protein;
    protein.defrost(); 
    protein.atoms().freeze(new AtomList.BackboneFilter());    
    distanceMatrix = new DistanceMatrix(protein.atoms(), 5.5, 2.0,4); 

    double[] tmprot;

    double[][] pp = new double[5000][3];
    double[][] targets = new double[5000][4];
    for (int c=0 ; c<targets.length ; c++)
       targets[c][0] = targets[c][1] = targets[c][2] = targets[c][3] = -999;
    phipsi(pp);

    // putting Rot1 into targets
    for (int kk=0 ; kk<pp.length ; kk++) {
    	if ((pp[kk][2]>-900) && (pp[kk][2]!=0) && (pp[kk][2]!=5)) {
    	    tmprot = lib.getRotamer((int) pp[kk][2] , pp[kk][0] , pp[kk][1] , 0);
    	    for (int dd=0 ; dd<lib.getChiMax((int) pp[kk][2]) ; dd++)
    	      targets[kk][dd] = tmprot[dd];
    	}
    }
    for (int kk=0 ; kk<protein.residues().size() ; kk++) 
		if ((protein.residues().residueAt(kk).type > -1) &
		    (protein.residues().residueAt(kk).type != UNK) &
		    residuesToModel.contains(protein.residues().residueAt(kk))){
		    ResidueBuilder.build(protein.residues().residueAt(kk), 
					 protein.residues().residueAt(kk).type, 
					 targets[kk]);
		}
    }


    protected void phipsi(double[][] pp) {
    for (int c=0 ; c<pp.length ; c++) {
       pp[c][0] = pp[c][1] = -60*Math.PI/180.0;
       pp[c][2]= -999;
    }
    TorsionList phiList = (TorsionList) TorsionList.createTorsionList(protein,
                distanceMatrix).filter(new TorsionList.FilterPhi(),
                                       new TorsionList());
    TorsionList psiList = (TorsionList) TorsionList.createTorsionList(protein,
                distanceMatrix).filter(new TorsionList.FilterPsi(),
                                       new TorsionList());
    for(int i=0 ; i<phiList.size() ; i++) {
        pp[((Torsion) phiList.elementAt(i)).getTorsionResNum()][0] =
           ((Torsion) phiList.elementAt(i)).torsion();
        pp[((Torsion) phiList.elementAt(i)).getTorsionResNum()][2] = Residue.type(
           ((Torsion) phiList.elementAt(i)).getTorsionResName());
    }
    for(int i=0 ; i<psiList.size() ; i++) {
        pp[((Torsion) psiList.elementAt(i)).getTorsionResNum()][1] =
           ((Torsion) psiList.elementAt(i)).torsion();
        pp[((Torsion) psiList.elementAt(i)).getTorsionResNum()][2] = Residue.type(
           ((Torsion) psiList.elementAt(i)).getTorsionResName());
    }
    }

} // of class

