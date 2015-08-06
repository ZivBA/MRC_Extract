package meshi.applications.homologyModelingNir;

import meshi.geometry.ResidueBuilder;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.residuesExtendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.Residues;

public class TrivialHomologyModeling implements Residues {
	
	
	/**
	 * In the returned array index n is a pointer to a single cell array with the  residue number 
	 * in the template for the query residue n. If there is no matching template residue the index
	 * points to null.
	 */
	public static int[][] alignmentMapping(int firstResInTemplate,
			String queryAlignment, String templateAlignment) {
		/* QNS1 = Query Not Starting in 1 */
		
		int queryLength = 0;
		for (int c=0 ; c<queryAlignment.length() ; c++)
			if (queryAlignment.charAt(c)!='-')
				queryLength++;
		int[][] result = new int[queryLength+1][]; /*QNS1*/
		
		int resNumCounterQ=0;
		int resNumCounterT=0;
		
		// Building the protein instance of the trivial homology
		resNumCounterQ = 0; /*QNS1*/
		resNumCounterT = firstResInTemplate-1;
		for (int position=0; position<queryAlignment.length() ; position++) {
			if (queryAlignment.charAt(position)!='-') 
				resNumCounterQ++;
			if (templateAlignment.charAt(position)!='-') 
				resNumCounterT++;		
			if (queryAlignment.charAt(position)!='-') {
				if (templateAlignment.charAt(position)!='-') {
					int[] tmp = new int[1];
					tmp[0] = resNumCounterT;
					result[resNumCounterQ] = tmp;
				}
			}
		}
		
		return result;
	}
	
	public static Protein trivialHomology(Protein template, int firstResInTemplate,
			String queryAlignment, String templateAlignment,  double[][] ppTemplate, DunbrackLib lib) {
		
		/* QNS1 = Query Not Starting in 1 */
		
		Protein query = null;
		int resNumCounterQ=0;
		int resNumCounterT=0;
		boolean[] completedResidue;
		double[][] ppQuery = new double[/*QNS1*/queryAlignment.length()+1][];

//		System.out.println("Alignment in template starts at: " + firstResInTemplate);
//		System.out.println("QUERY " + queryAlignment);
//		System.out.println("TEMPL " + templateAlignment);
		completedResidue = new boolean[queryAlignment.length()];
		for (int c=0; c<queryAlignment.length() ; c++)
			completedResidue[c] = true;

		
		// Building the protein instance of the trivial homology
		AtomList CAslist = new AtomList();
		resNumCounterQ = 0; /*QNS1*/
		resNumCounterT = firstResInTemplate-1;
		for (int position=0; position<queryAlignment.length() ; position++) {
			if (queryAlignment.charAt(position)!='-') 
				resNumCounterQ++;
			if (templateAlignment.charAt(position)!='-') 
				resNumCounterT++;		
			if (queryAlignment.charAt(position)!='-') {
				if ((templateAlignment.charAt(position)!='-') && (template.residue(resNumCounterT)!= null) && 
						(template.residue(resNumCounterT).ca()!= null)) {
					CAslist.add(new Atom(0,0,0,"CA",
							Residue.one2three(queryAlignment.charAt(position)),
							resNumCounterQ,-1));
					ppQuery[resNumCounterQ] = ppTemplate[resNumCounterT];
				}
			}
		}
		query = new Protein(CAslist, new ResidueExtendedAtoms(ADD_ATOMS));
				 
		// Assigning the coordinates to the trivial homology
		resNumCounterQ = 0; /*QNS1*/
		resNumCounterT = firstResInTemplate-1;
		for (int position=0; position<queryAlignment.length() ; position++)  {	
			if (queryAlignment.charAt(position)!='-') 
				resNumCounterQ++;
			if (templateAlignment.charAt(position)!='-') 
				resNumCounterT++;		
			if ((queryAlignment.charAt(position)!='-') && (templateAlignment.charAt(position)!='-') &&
				(template.residue(resNumCounterT)!= null) && (template.residue(resNumCounterT).ca()!= null)) {
				for (int atomCounter=0 ; atomCounter<query.residue(resNumCounterQ).atoms().size() ; atomCounter++) {
//					System.out.println(resNumCounterQ + " " + queryAlignment.charAt(position) + " " + resNumCounterT + " " + templateAlignment.charAt(position));
					Atom atomInQuery = query.residue(resNumCounterQ).atoms().atomAt(atomCounter);
					Atom atomInTemplate = template.residue(resNumCounterT).atoms().getAtom(atomInQuery.name());
					if (atomInTemplate!=null) {
						atomInQuery.setXYZ(atomInTemplate.x(), atomInTemplate.y(), atomInTemplate.z());
						atomInQuery.freeze();
					}
					else {
						if (!atomInQuery.isHydrogen) // We are missing a heavy atom in the template
							completedResidue[position] = false;
						atomInQuery.defrost();
					}
				}
			}
		}
		
		// Brining new atoms close to their bonded atoms
		for (int atomCounter=0 ; atomCounter<query.atoms().size() ; atomCounter++) {
			Atom atom = query.atoms().atomAt(atomCounter);
			if (!atom.frozen())
				atom.setXYZ(atom.bonded().atomAt(0).x()+Math.random(), atom.bonded().atomAt(0).y()+Math.random(), atom.bonded().atomAt(0).z()+Math.random());
		}

/*
 *  I don't need that now:
 *  ----------------------
		// First Minimization - This will mainly set the hydrogens right.
		EnergyCreator[] energyCreators = {  
				new BondCreator(),
				new AngleCreator(),
				new PlaneCreator(),
				new OutOfPlaneCreator()
		};
		DistanceMatrix distanceMatrix = new DistanceMatrix(query.atoms(), 5.5,  2.0,  4);  
		TotalEnergy energy = new TotalEnergy(query, distanceMatrix, energyCreators, commands);
		Minimizer minimizer = new LBFGS(energy, 0.001 , 100000 , 100 );  
		System.out.println(minimizer.minimize());
		
		// Freezing only residues whose sidechains are determined from the template
		resNumCounterQ = 0;
		for (int position=0; position<queryAlignment.length() ; position++) 
			if (queryAlignment.charAt(position)!='-') {
				resNumCounterQ++;
				if ((templateAlignment.charAt(position)==queryAlignment.charAt(position)) &&
						completedResidue[position]) { 
					//System.out.println("Freezing: " + resNumCounterQ + " " + templateAlignment.charAt(position));
					query.residue(resNumCounterQ).atoms().freeze();
				}
				else {
					if (query.residue(resNumCounterQ)!=null) {
						//System.out.println("Defrost: " + resNumCounterQ + " " + templateAlignment.charAt(position));
						query.residue(resNumCounterQ).atoms().defrost();
					}
				}
			}
		
		// Running SCMOD
		SCMOD.scmod(commands, new DunbrackLib(commands,1.0,50), query, 2);
*/		
		
		// Putting query in ROT1		
		for (int res = 0 ; res<query.residues().size() ; res++)
			if (query.residues().residueAt(res).ca()!=null)
				if ( (query.residues().residueAt(res).type!=GLY) && (query.residues().residueAt(res).type!=ALA) ) {
					int inRes = query.residues().residueAt(res).number;
					int resType = query.residues().residueAt(res).type;
					ResidueBuilder.build(query.residue(inRes),  resType, 
							lib.getRotamer(resType, ppQuery[inRes][0], ppQuery[inRes][1], 0));
					}
				else {
					int inRes = query.residues().residueAt(res).number;
					int resType = query.residues().residueAt(res).type;
					ResidueBuilder.build(query.residue(inRes),  resType, null);				
				}
		
		query.defrost();
		return query;
	}

}
