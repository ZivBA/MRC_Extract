package meshi.applications.prediction;

import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.angle.AngleCreator;
import meshi.energy.bond.BondCreator;
import meshi.energy.linearRG.LinearRgCreator;
import meshi.energy.outOfPlane.OutOfPlaneCreator;
import meshi.energy.plane.PlaneCreator;
import meshi.energy.softExcludedVol.SoftExcludedVolCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LBFGS;
import meshi.optimizers.Minimizer;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.sequences.Sequence;
import meshi.sequences.SequenceAlignment;
import meshi.sequences.SequenceAlignmentColumn;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.rotamericTools.RotamericTools;


public class CompleteClass implements Residues, AtomTypes , KeyWords { 

 
 	/**
 	* Completing the all the missing residues in the protein.
 	**/
    public static Protein complete(CommandList commands , Protein holed) throws Exception {
    	Sequence sequence1 = new PredictionResidueSequence(commands, TARGET_FILE_PATH);
    	Sequence sequence2 = holed.sequence();
    	SequenceAlignment alignment = SequenceAlignment.identityAlignment(sequence1, sequence2);
    	return complete(commands , holed , 1 , ((SequenceAlignmentColumn) alignment.columnAt(alignment.size()-1)).cell0().number);     	
    }


 	/**
 	* Completing the missing residues in the interval [firstToComplete,lastToComplete] of the protein.
 	* Followed is a short minimization where the gap edge residues are also defrosted.
 	* The protein is put into Rot1 and returned.
 	*
 	* All the hydrogens are also added, and minimized to the right place.
 	**/
    public static Protein complete(CommandList commands , Protein holed , int firstToComplete, int lastToComplete) throws Exception {
    	Sequence sequence1 = new PredictionResidueSequence(commands, TARGET_FILE_PATH);
    	Sequence sequence2 = holed.sequence();
    	SequenceAlignment alignment = SequenceAlignment.identityAlignment(sequence1, sequence2);
    	System.out.println("The holed protein aligned with the full sequence:"+"\n"+alignment);
    	AtomList newList = new AtomList();
    	int prev = findNonDummy(holed,firstToComplete-1);
    	int next = findNonDummy(holed,prev);
    	double bx=1e10,by=1e10,bz=1e10,ex=1e10,ey=1e10,ez=1e10;
    	int[] added = new int[((SequenceAlignmentColumn) alignment.columnAt(alignment.size()-1)).cell0().number+1];
    	for (int i=1 ; i<added.length ; i++) {
    		if (((SequenceAlignmentColumn) alignment.getColumn(0,i)).getChar(0) != 'X') { // Skipping the dummy residues in the begining
    			if  ((i<firstToComplete) || (i>lastToComplete)) {  // Not in the completeing interval
    				if (holed.residue(i) != null)
	    				for (int j=0 ; j<holed.residue(i).atoms().size() ; j++)
	    					newList.add(new Atom(holed.residue(i).atoms().atomAt(j)));
    			}
    			else { //in the completeing interval 
	    			if (!((SequenceAlignmentColumn) alignment.getColumn(0,i)).cell1().gap()) { // no need to complete
    					for (int j=0 ; j<holed.residue(i).atoms().size() ; j++)
    						newList.add(new Atom(holed.residue(i).atoms().atomAt(j)));
    					if (i>=next) { // need to update prev and next
    						int tmpind = findNonDummy(holed,next);
    						if (tmpind>0) {
    							prev = next;
    							next = tmpind;
    						}
    					}
    				}
    				else { // Have to complete
						bx = holed.residue(prev).ca().x();
						by = holed.residue(prev).ca().y();
						bz = holed.residue(prev).ca().z();			
						ex = holed.residue(next).ca().x();
						ey = holed.residue(next).ca().y();
						ez = holed.residue(next).ca().z();
						newList.add(new Atom(bx+(ex-bx)*(i-prev)*1.0/(next-prev),
											 by+(ey-by)*(i-prev)*1.0/(next-prev),
											 bz+(ez-bz)*(i-prev)*1.0/(next-prev),
											 "CA", 
											 Residue.one2three("" + ((SequenceAlignmentColumn) alignment.getColumn(0,i)).getChar(0)),
											 i,
											 -1));
						System.out.println("Added: " + newList.atomAt(newList.size()-1));
						added[i] = 1;
					}
				}
			}
		}
		
		// Filling the rest of atoms, and freezeing the correct part		
		Protein pr = new Protein(newList,new ResidueExtendedAtoms(ADD_ATOMS));
		pr.defrost();
		pr.freeze(); //pr.atoms().filter(new AtomList.NonHydrogen()).freeze();
		for (int i=1 ; i<=pr.residues().residueAt(pr.residues().size()-1).number ; i++)
			if (pr.residue(i)!=null) {
				if (added[i]>0) {
					if (pr.residue(i-1)!=null)
						pr.residue(i-1).atoms().defrost();
					if (pr.residue(i)!=null)
						pr.residue(i).atoms().defrost();
					if (pr.residue(i+1)!=null)
						pr.residue(i+1).atoms().defrost();
				}
			}

		System.out.println("************* Minimizing *********");
   		DistanceMatrix dm = new DistanceMatrix(pr.atoms(),  5.5, 2.0, 3);  // Note that it is 3 bonds not 4
		EnergyCreator[] energyCreators = {  
	    	new BondCreator(),
	    	new AngleCreator(),
	    	new PlaneCreator(),
	    	new OutOfPlaneCreator(),
	   	 	new SoftExcludedVolCreator(2.0,2),
	    	new LinearRgCreator(10.0)
		};
		TotalEnergy energy = new TotalEnergy(pr, dm, energyCreators, commands);
		Minimizer minimizer = new LBFGS(energy, 0.5 , 20000 , 100);  
		System.out.println(minimizer.minimize());

		System.out.println("************* Into Rot 1 *********");
   		pr.defrost();
   		dm = new DistanceMatrix(pr.atoms(),  5.5, 2.0, 4);
   		DunbrackLib lib = new DunbrackLib(commands , 1.0 , 2);
   		RotamericTools.putIntoRot1(pr, dm, lib);
   		return pr;
   	}

    protected static int findNonDummy(Protein holed , int begin) {
    	if (begin>=holed.residues().residueAt(holed.residues().size()-1).number)
    		return -1;
    	for (int i=begin+1 ; i<=holed.residues().residueAt(holed.residues().size()-1).number ; i++)
    		if ((holed.residue(i) != null) && !holed.residue(i).dummy())  
    			return i;
    	return -1;
    }

}
