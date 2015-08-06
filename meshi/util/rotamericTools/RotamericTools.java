package meshi.util.rotamericTools;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.ResidueBuilder;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionList;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.parameters.AtomTypes;
import meshi.parameters.Residues;
import meshi.util.KeyWords;
import meshi.util.overlap.Overlap;


public class RotamericTools implements Residues, AtomTypes , KeyWords, Rot1Arrays { 


	// This is a simplified version of the next overload of the method. It assumes that pp already exists.
	// Also the distance matrix is not discarded.
	public static void putIntoRot1(Protein prot, double[][] pp, DunbrackLib lib) { 

		for (int kk=0 ; kk<pp.length ; kk++) 
			if ((pp[kk]!=null) && (pp[kk][2]>-1) && (prot.residue(kk).atoms().defrostedAtoms().size()>0)) {
				if (pp[kk][2] ==ALA)
					ResidueBuilder.build(prot.residue(kk),0,null);
				else if (pp[kk][2]!=GLY) { 
					ResidueBuilder.build(prot.residue(kk),prot.residue(kk).type,
							lib.getRotamer((int) pp[kk][2] , pp[kk][0] , pp[kk][1] , 0));
				}
			}    
	}  // of putIntoRot1

	// This method puts all the residues of prot, that have at least ONE NON-FROZEN atom, into their rot1 position.
	// It also returns an Nx4 double array where row n is the phi,psi,type,{rotamer 1 prob} of residue n.
	// Note that the distance matrix is discarded after this method is run, as it is the best way to ensure,
	// that it is updated properly after this sidechain change.
	public static double[][] putIntoRot1(Protein prot, DistanceMatrix dm, DunbrackLib lib) { 

		double[][] pp = phipsi(prot,dm);
		double[] tmprot;

		int maxRes = -1;
		for (int kk=0 ; kk<prot.atoms().size() ; kk++) 
			if (prot.atoms().atomAt(kk).residueNumber() > maxRes)
				maxRes = prot.atoms().atomAt(kk).residueNumber();
		double[][] largepp = new double[maxRes+1][];
		for (int kk=0 ; kk<largepp.length ; kk++) 
			if ((prot.residue(kk)!=null) && (prot.residue(kk).type>-1))
				if ((kk<pp.length) && (pp[kk]!=null))
					largepp[kk] = pp[kk];
				else {
					largepp[kk] = new double[3];
					largepp[kk][0] = -60*Math.PI/180.0;
					largepp[kk][1] = -40*Math.PI/180.0;
					largepp[kk][2] = prot.residue(kk).type;
				}
		pp = largepp;		

		// putting Rot1 into all instances of this residue type
		for (int kk=0 ; kk<pp.length ; kk++) {
			if ((pp[kk]!=null) && (pp[kk][2]>-1) && (prot.residue(kk).atoms().defrostedAtoms().size()>0)) {
				if (pp[kk][2] ==ALA)
					ResidueBuilder.build(prot.residue(kk),0,null);
				else if (pp[kk][2]!=GLY) {
					tmprot = lib.getRotamer((int) pp[kk][2] , pp[kk][0] , pp[kk][1] , 0);
					ResidueBuilder.build(prot.residue(kk),prot.residue(kk).type,tmprot);
				}
			}
		}

		double[][] tmppp = new double[pp.length][];
		for (int kk=0 ; kk<pp.length ; kk++) {
			if ((pp[kk]!=null) && (pp[kk][2]>ALA) && (pp[kk][2]!=GLY) &&
					(prot.residue(kk).atoms().defrostedAtoms().size()>0)) {
				tmppp[kk] = new double[4];
				tmppp[kk][0] = pp[kk][0];
				tmppp[kk][1] = pp[kk][1];
				tmppp[kk][2] = pp[kk][2];
				tmppp[kk][3] = lib.getRotamerProb((int) pp[kk][2] , pp[kk][0] , pp[kk][1] , 0); 	
			}
			else if ((pp[kk]!=null) && ((pp[kk][2]==ALA) || (pp[kk][2]==GLY)) &&
					(prot.residue(kk).atoms().defrostedAtoms().size()>0)) {
				tmppp[kk] = new double[3];
				tmppp[kk][0] = pp[kk][0];
				tmppp[kk][1] = pp[kk][1];
				tmppp[kk][2] = pp[kk][2];
			}
		}
		// Discarding the distance matrix.
//		dm.terminator.kill("The distance matrix is not updated correctly in PutIntoRot1. It has to be created anew.");
		return tmppp;    
	}  // of putIntoRot1

	
	/** This method puts all the residues of prot, that have at least ONE NON-FROZEN atom, into their 
	 * most similar rotamer. It also returns an Nx8 double array where row n is the:
	 * phi,psi,type,{rotamer prob}, {chi1}, {chi2}, {chi3}, {chi4}, {dev1}, {dev2}, {dev3}, {dev4} of residue n.
	 **/
	public static double[][] putIntoNearestRot(Protein prot, DistanceMatrix dm, DunbrackLib lib, boolean changeRot) { 

		double[][] pp = phipsi(prot,dm);
		double[] tmprot;
		double[] tmprotdev;

		// Building the array to return
		int maxRes = -1;
		for (int kk=0 ; kk<prot.atoms().size() ; kk++) 
			if (prot.atoms().atomAt(kk).residueNumber() > maxRes)
				maxRes = prot.atoms().atomAt(kk).residueNumber();
		double[][] largepp = new double[maxRes+1][];
		for (int kk=0 ; kk<largepp.length ; kk++) 
			if ((prot.residue(kk)!=null) && (prot.residue(kk).type>-1))
				if ((kk<pp.length) && (pp[kk]!=null)) {
					largepp[kk] = new double[12];
					largepp[kk][0] = pp[kk][0];
					largepp[kk][1] = pp[kk][1];
					largepp[kk][2] = prot.residue(kk).type;
				}
				else {
					largepp[kk] = new double[12];
					largepp[kk][0] = -60*Math.PI/180.0;
					largepp[kk][1] = -40*Math.PI/180.0;
					largepp[kk][2] = prot.residue(kk).type;
				}
		pp = largepp;		

		double[][] chis = getSidechainTorsions(prot);
		for (int kk=0 ; kk<pp.length ; kk++) {
			Residue residue = prot.residue(kk);
			if ((pp[kk]!=null) && (pp[kk][2]>-1) && (residue.atoms().defrostedAtoms().size()>0)) {
				if (pp[kk][2] ==ALA)
					ResidueBuilder.build(residue,0,null);
				else if (pp[kk][2]!=GLY) {
					int closestRot = nearestRotAngles(lib, residue.type, pp[kk][0], pp[kk][1], chis[kk], 45.0);
					if (closestRot==-1) {
						System.out.println("Change Rot in: " + residue);
						closestRot = nearestRotReturningRotamerInd(lib, residue, pp[kk][0] , pp[kk][1]);
					}
//					System.out.println(kk + " " + Residue.nameOneLetter(residue.type) + " " + closestRot);
					tmprot = lib.getRotamer(residue.type , pp[kk][0] , pp[kk][1] , closestRot);
					tmprotdev = lib.getRotamerDev(residue.type , pp[kk][0] , pp[kk][1] , closestRot);
					if (changeRot) {
						ResidueBuilder.build(prot.residue(kk),prot.residue(kk).type,tmprot);
					}
					pp[kk][3] = lib.getRotamerProb(residue.type , pp[kk][0] , pp[kk][1] , closestRot);
					for (int chiCounter=0 ; chiCounter<tmprot.length ; chiCounter++) {
						pp[kk][4+chiCounter] = tmprot[chiCounter];
						pp[kk][8+chiCounter] = tmprotdev[chiCounter];
					}
				}
			}
		}
		return pp;    
	} 


	// Returns the an Nx3 double array where row n is the phi,psi,type of residue n.
	public static double[][] phipsi(Protein prot, DistanceMatrix dm) {
		int maxResNum = -999;
		int resNum;
		double[][] pp;

		TorsionList phiList = (TorsionList) TorsionList.createTorsionList(prot,
				dm).filter(new TorsionList.FilterPhi(),
						new TorsionList());
		TorsionList psiList = (TorsionList) TorsionList.createTorsionList(prot,
				dm).filter(new TorsionList.FilterPsi(),
						new TorsionList());

		for(int i=0 ; i<phiList.size() ; i++) 
			if (((Torsion) phiList.elementAt(i)).getTorsionResNum()>maxResNum) 
				maxResNum = ((Torsion) phiList.elementAt(i)).getTorsionResNum();

		for(int i=0 ; i<psiList.size() ; i++) 
			if (((Torsion) psiList.elementAt(i)).getTorsionResNum()>maxResNum) 
				maxResNum = ((Torsion) psiList.elementAt(i)).getTorsionResNum();

		pp = new double[maxResNum+1][];

		for(int i=0 ; i<phiList.size() ; i++) {
			resNum = ((Torsion) phiList.elementAt(i)).getTorsionResNum();
			if (pp[resNum] == null) {
				pp[resNum] = new double[3];
				pp[resNum][0] = -60*Math.PI/180.0;
				pp[resNum][1] = -40*Math.PI/180.0;
				pp[resNum][2] = -1;
			}
			pp[resNum][0] = ((Torsion) phiList.elementAt(i)).torsion();
			pp[resNum][2] = Residue.type(((Torsion) phiList.elementAt(i)).getTorsionResName());
		}
		for(int i=0 ; i<psiList.size() ; i++) {
			resNum = ((Torsion) psiList.elementAt(i)).getTorsionResNum();
			if (pp[resNum] == null) {
				pp[resNum] = new double[3];
				pp[resNum][0] = -60*Math.PI/180.0;
				pp[resNum][1] = -40*Math.PI/180.0;
				pp[resNum][2] = -1;
			}
			pp[resNum][1] = ((Torsion) psiList.elementAt(i)).torsion();
			pp[resNum][2] = Residue.type(((Torsion) psiList.elementAt(i)).getTorsionResName());
		}

		return pp;
	} // of phipsi



	// This method finds the most similar (in term of rms) rotamer of residue 'res'.
	// If it returns -1 then an error occured, probably some atoms are missing in the side chain. 
	// Otherwise it returns the rotamer index found. 
	public static int nearestRotReturningRotamerInd(DunbrackLib lib, Residue res, double phi, double psi) { 
		AtomList al1 = res.atoms();
		AtomList al2 = al1.duplicate();
		int restype = res.type;
		int numOfRotamers = lib.getRotamerNum(restype , phi , psi);
		int bestrot = -1;
		double[] tmprot;
		double rms,minrms = Double.MAX_VALUE;
		boolean cont = true;
		int ind=0;
		for ( ; cont && (ind<numOfRotamers) ; ind++) {
			tmprot = lib.getRotamer(restype , phi , psi , ind);
			try {
				ResidueBuilder.build(res, restype, tmprot);
			}
			catch (Exception e) {
				cont = false;
			}
			rms = calcRMS(al1,al2);
			if (rms<0.0) 
				cont = false;
			else if (rms<minrms) {
				minrms=rms;
				bestrot = ind;
			}
			if ((restype==ASP) || (restype==GLU) || (restype==PHE) || (restype==ARG) || (restype==TYR)) {
				try {
					flipTerm(res);
				}
				catch (Exception e) {
					cont = false;
				}
				if ((restype==ASP) || (restype==PHE) || (restype==TYR)) {
					tmprot[1] += Math.PI;
					if (tmprot[1] > Math.PI)
						tmprot[1] -= 2*Math.PI;
				}
				if (restype==GLU) {
					tmprot[2] += Math.PI;
					if (tmprot[2] > Math.PI)
						tmprot[2] -= 2*Math.PI;
				}
				rms = calcRMS(al1,al2);
				if (rms<0.0) 
					cont = false;
				else if (rms<minrms) {
					minrms=rms;
					bestrot = ind;
				}                  
				try {
					flipTerm(res);
				}
				catch (Exception e) {
					cont = false;
				}
			}
		}
		// returning the positions in residue to the original
		for (int c=0 ; c<al2.size() ; c++) {
			if (!al1.atomAt(c).name.equals(al2.atomAt(c).name)) {
				throw new RuntimeException("The lists are not ordered the same. This should not have happened!!");
			}
			al1.atomAt(c).setXYZ(al2.atomAt(c).x(),al2.atomAt(c).y(),al2.atomAt(c).z());
		}
		if (cont)
			return bestrot;
		else
			return -1;
	}
	
	
	// This method finds the most similar (in term of rms) rotamer of residue 'res'.
	// If it returns an array then a similar rotamer was found in the Dunbrack library, and it is returned. 
	// If it returns null then an error occured, probably some atoms are missing in the side chain. 
	public static double[] nearestRot(DunbrackLib lib, Residue res, double phi, double psi) {
		int ind = nearestRotReturningRotamerInd(lib, res, phi, psi);
		if (ind == -1)
			return null;
		else 
			return lib.getRotamer(res.type , phi , psi , ind);
	}


	// al1,al2 - are two atom lists from two instances a certain residue. 
	// This method overlap the backbone atoms (N,C,CA,O) of the two residues and returns the RMS of their
	// (non-hydrogen) sidechain atoms (including the CB atom).
	// If the two lists do not contain exactly the same atoms (in term of names not instances) then -1 is returned.
	// The order in the lists is irrelevent.
	public static double calcRMS(AtomList al1, AtomList al2) {
		al1 = al1.filter(new AtomList.NonHydrogen()).noOXTFilter();
		al2 = al2.filter(new AtomList.NonHydrogen()).noOXTFilter();
		if (al1.size()!=al2.size())
			return -1;
		boolean found;
		double rms = 0.0;
		int NUM_BACKBONE_ATOMS = 4;
		int backboneCount = 0;
		int anum = 0;
		int arCount = NUM_BACKBONE_ATOMS;
		double[][] co = new double[3][al1.size()];
		double[][] co2 = new double[3][al1.size()];
		int[] aux = new int[NUM_BACKBONE_ATOMS];
		for (int c=0; c<NUM_BACKBONE_ATOMS ; c++)
			aux[c] = c;
		for (int c1=0 ; c1<al1.size() ; c1++) {
			found = false;
			for (int c2=0 ; c2<al2.size() ; c2++)
				if (al1.atomAt(c1).name.equals(al2.atomAt(c2).name)) {
					found = true;
					if (al1.atomAt(c1).name.equals("N")) {
						co[0][0] = al1.atomAt(c1).x();
						co[1][0] = al1.atomAt(c1).y();
						co[2][0] = al1.atomAt(c1).z();
						co2[0][0] = al2.atomAt(c2).x();
						co2[1][0] = al2.atomAt(c2).y();
						co2[2][0] = al2.atomAt(c2).z();
						backboneCount++;
					}
					else if (al1.atomAt(c1).name.equals("CA")) {
						co[0][1] = al1.atomAt(c1).x();
						co[1][1] = al1.atomAt(c1).y();
						co[2][1] = al1.atomAt(c1).z();
						co2[0][1] = al2.atomAt(c2).x();
						co2[1][1] = al2.atomAt(c2).y();
						co2[2][1] = al2.atomAt(c2).z();
						backboneCount++;
					}
					else if (al1.atomAt(c1).name.equals("C")) {
						co[0][2] = al1.atomAt(c1).x();
						co[1][2] = al1.atomAt(c1).y();
						co[2][2] = al1.atomAt(c1).z();
						co2[0][2] = al2.atomAt(c2).x();
						co2[1][2] = al2.atomAt(c2).y();
						co2[2][2] = al2.atomAt(c2).z();
						backboneCount++;
					}
					else if (al1.atomAt(c1).name.equals("O")) {
						co[0][3] = al1.atomAt(c1).x();
						co[1][3] = al1.atomAt(c1).y();
						co[2][3] = al1.atomAt(c1).z();
						co2[0][3] = al2.atomAt(c2).x();
						co2[1][3] = al2.atomAt(c2).y();
						co2[2][3] = al2.atomAt(c2).z();
						backboneCount++;
					}
					/*else if (al1.atomAt(c1).name.equals("CB")) {  You should revive this section and update NUM_BACKBONE_ATOMS by one, if you want to use the CB in the alignment
    	      	 	anum++;
    	      	 	co[0][4] = al1.atomAt(c1).x();
    	      	 	co[1][4] = al1.atomAt(c1).y();
    	      	 	co[2][4] = al1.atomAt(c1).z();
    	      	 	co2[0][4] = al2.atomAt(c2).x();
    	      	 	co2[1][4] = al2.atomAt(c2).y();
    	      	 	co2[2][4] = al2.atomAt(c2).z();
    	      	 	backboneCount++;
    	      	 }  */
					else {
						if (arCount>=co[0].length)
							return -1;
						anum++;   
						co[0][arCount] = al1.atomAt(c1).x();
						co[1][arCount] = al1.atomAt(c1).y();
						co[2][arCount] = al1.atomAt(c1).z();
						co2[0][arCount] = al2.atomAt(c2).x();
						co2[1][arCount] = al2.atomAt(c2).y();
						co2[2][arCount] = al2.atomAt(c2).z();
						arCount++;    	      	 	
					}
				}
			if (!found)
				return -1;
		}
		if (backboneCount<NUM_BACKBONE_ATOMS) {
			return -1;
		}

		Overlap.rmsPartial(co,co2,aux);
		for (int c1=NUM_BACKBONE_ATOMS-1; c1<co[0].length ; c1++) 
			rms = rms + (co[0][c1] - co2[0][c1])*(co[0][c1] - co2[0][c1]) + 
			(co[1][c1] - co2[1][c1])*(co[1][c1] - co2[1][c1]) + 
			(co[2][c1] - co2[2][c1])*(co[2][c1] - co2[2][c1]);

		rms = rms/anum;
		return rms;
	}


	// This method returns a Nx4 array with the sidechain torsion values of the protein.
	// N is the residue number of the last res in the protein. In each row the values are {chi1,chi2,chi3,chi4}
	// This way [6][2] is the chi3 value of residue 6.
	public static double[][] getSidechainTorsions(Protein protein) {
		DistanceMatrix	distanceMatrix = new DistanceMatrix(protein.atoms(), 5.5, 2.0,4);
		double[][] targets = new double[protein.residues().residueAt(protein.residues().size()-1).number+1][4];
		try {
			distanceMatrix.update(1); 
		}
		catch (Exception e) {
			throw new RuntimeException(e);
		}
		TorsionList modelChi1 = (TorsionList) TorsionList.createTorsionList(protein,distanceMatrix).filter(new TorsionList.FilterChi1() , new TorsionList());
		TorsionList modelChi2 = (TorsionList) TorsionList.createTorsionList(protein,distanceMatrix).filter(new TorsionList.FilterChi2() , new TorsionList());
		TorsionList modelChi3 = (TorsionList) TorsionList.createTorsionList(protein,distanceMatrix).filter(new TorsionList.FilterChi3() , new TorsionList());
		TorsionList modelChi4 = (TorsionList) TorsionList.createTorsionList(protein,distanceMatrix).filter(new TorsionList.FilterChi4() , new TorsionList());
		for (int cc=0 ; cc<modelChi1.size() ; cc++) 
			targets[modelChi1.torsionAt(cc).getTorsionResNum()][0] = modelChi1.torsionAt(cc).torsion();
		for (int cc=0 ; cc<modelChi2.size() ; cc++) 
			targets[modelChi2.torsionAt(cc).getTorsionResNum()][1] = modelChi2.torsionAt(cc).torsion();
		for (int cc=0 ; cc<modelChi3.size() ; cc++) 
			targets[modelChi3.torsionAt(cc).getTorsionResNum()][2] = modelChi3.torsionAt(cc).torsion();
		for (int cc=0 ; cc<modelChi4.size() ; cc++) 
			targets[modelChi4.torsionAt(cc).getTorsionResNum()][3] = modelChi4.torsionAt(cc).torsion();
		return targets;
	}

	// This method returns the number of the highest probability rotamer in the library whose angles are the first 
	// to be all within "limit" of the residue's. "limit" should not be more than 30-40 degrees.
	// If no rotamer was found to be within the "limit"s - (-1) is returned.   
	// Limit should be given in DEGREES.
	// One of the inputs to the method is the sidechain angles. The method "getSidechainTorsions" in this class can get them for you. 
	// Any "restype" not in [1..4,6..19] (Ala,Gly,Unk etc.) will return (-1) too.
	public static int nearestRotAngles(DunbrackLib lib, int restype, double phi, double psi, 
			double[] sidechainAng, double limit) {
		if ((restype<CYS) || (restype==GLY) || (restype>TYR))
			return -1;
		limit = limit*Math.PI/180.0;
		int numOfRotamers = lib.getRotamerNum(restype , phi , psi);
		int maxChi = lib.getChiMax(restype);
		double[] tmprot;
		double tmpVal;
		boolean inlimit = false;
		for (int ind=0 ; ind<numOfRotamers ; ind++) {
			tmprot = lib.getRotamer(restype , phi , psi , ind);
			inlimit = true;
			for (int d=0 ; d<maxChi ; d++) {
				tmpVal = Math.abs(sidechainAng[d]-tmprot[d]);
				if (tmpVal>Math.PI)
					tmpVal = 2*Math.PI - tmpVal;
				if (tmpVal>limit)
					inlimit=false;
			}
			if (inlimit)
				return ind;
			// Flipping symetrical plane in the side chain.
			if ((restype==ASP) || (restype==PHE) || (restype==TYR)) {
				tmprot[1] += Math.PI;
				if (tmprot[1]>Math.PI)
					tmprot[1] -= 2*Math.PI;
			}
			if (restype==GLU) {
				tmprot[2] += Math.PI;
				if (tmprot[2]>Math.PI)
					tmprot[2] -= 2*Math.PI;
			}
			inlimit = true;
			for (int d=0 ; d<maxChi ; d++) {
				tmpVal = Math.abs(sidechainAng[d]-tmprot[d]);
				if (tmpVal>Math.PI)
					tmpVal = 2*Math.PI - tmpVal;
				if (tmpVal>limit)
					inlimit=false;
			}
			if (inlimit)
				return ind;
		}
		return -1;
	}

	// This methods assign every sidechain torsion a random angle.
	// All the atoms must be present in the protein. 
	public static void jumble(Protein prot){
		double[] tmprot={-999,-999,-999,-999};
		for (int kk=0; kk<prot.residues().size() ; kk++) {
			tmprot[0] = 2*Math.PI*(0.5-Math.random());
			tmprot[1] = 2*Math.PI*(0.5-Math.random());
			tmprot[2] = 2*Math.PI*(0.5-Math.random());
			tmprot[3] = 2*Math.PI*(0.5-Math.random());
			if ((prot.residue(kk).type>0) && (prot.residue(kk).type!=5) && (prot.residue(kk).type<20))
				ResidueBuilder.build(prot.residue(kk),prot.residue(kk).type,tmprot);
		}	
	} 


	// This method finds the most similar (in term of rms) rotamer of residue 'res', from a double array of 
	// angles , given in 'angles'. The first index in angles is the rotamer, and the second index runs on the chi's.
	// It returns the index of the most similar rotamer in 'angles'. 
	// If it returns -1 then an error occured, probably some atoms are missing in the side chain. 
	public static int nearestRotFromArray(Residue res, double[][] angles) { 
		AtomList al1 = res.atoms();
		AtomList al2 = al1.duplicate();
		int restype = res.type;
		int bestrot = -1;
		double[] tmprot;
		double rms,minrms = 999999999.9;
		boolean cont = true;
		for (int ind=0 ; cont && (ind<angles.length) ; ind++) {
			tmprot = new double[angles[ind].length];
			for (int tmp=0 ; tmp<tmprot.length ; tmp++) 
				tmprot[tmp] = angles[ind][tmp];
			try {
				ResidueBuilder.build(res, restype, tmprot);
			}
			catch (Exception e) {
				cont = false;
			}
			rms = calcRMS(al1,al2);
			if (rms<0.0) 
				cont = false;
			else if (rms<minrms) {
				minrms=rms;
				bestrot = ind;
			}
			if ((restype==ASP) || (restype==GLU) || (restype==PHE) || (restype==ARG) || (restype==TYR)) {
				try {
					flipTerm(res);
				}
				catch (Exception e) {
					cont = false;
				}
				rms = calcRMS(al1,al2);
				if (rms<0.0) 
					cont = false;
				else if (rms<minrms) {
					minrms=rms;
					bestrot = ind;
				}                  
				try {
					flipTerm(res);
				}
				catch (Exception e) {
					cont = false;
				}
			}
		}
		// returning the positions in residue to the original
		for (int c=0 ; c<al2.size() ; c++) {
			if (!al1.atomAt(c).name.equals(al2.atomAt(c).name)) {
				throw new RuntimeException("The lists are not ordered the same. This should not have happened!!");
			}
			al1.atomAt(c).setXYZ(al2.atomAt(c).x(),al2.atomAt(c).y(),al2.atomAt(c).z());
		}
		if (cont)
			return bestrot;
		else
			return -1;
	}




	public static double[] getMean(int type , double prob, int what) {
		double[] result = new double[2];
		int best=-999,secondbest=-999;
		for (int c=0 ; c<10 ; c++) 
			if (mean[type][c][what]>-0.999) {
				if (Math.abs(c/10.0+0.05 - prob) <= Math.abs(best/10.0+0.05 - prob)) {
					secondbest = best;
					best = c;
				}
				else if (Math.abs(c/10.0+0.05 - prob) <= Math.abs(secondbest/10.0+0.05 - prob)) {
					secondbest=c;    			
				}
			} 
		if (Math.abs(best/10.0+0.05 - prob) > 0.162)
			System.out.println("\n\nWarning!!! A prob of: " + prob + " that is too far from the closest:" + 
					(best/10.0+0.05) + " type: " + type + "\n\n");
		result[0] = mean[type][best][what] + 
		(mean[type][secondbest][what] - mean[type][best][what])/(secondbest/10.0 - best/10.0)*(prob-(best/10.0+0.05));
		result[1] = std[type][best][what] + 
		(std[type][secondbest][what] - std[type][best][what])/(secondbest/10.0 - best/10.0)*(prob-(best/10.0+0.05));

		return result;
	}


	// This method finds the rms between the side chain atoms of two residues of the same type.  
	// -1 is returned in case of an error.
	public static double calcRMS(Residue res1, Residue res2) { 
		int restype = res1.type;
		if ((restype<ALA) || (restype>TYR))
			return -1;
		if ((restype==ALA) || (restype==GLY))
			return 0.0;
		double rms = calcRMS(res1.atoms(),res2.atoms());
		double minrms = rms;

		if ((restype==ASP) || (restype==GLU) || (restype==PHE) || (restype==ARG) ||  (restype==TYR)) {
			flipTerm(res2);
			rms = calcRMS(res1.atoms(),res2.atoms());
			if ((rms<minrms) && (rms > -1))
				minrms=rms;
			flipTerm(res2);              
		}

		return minrms;
	}


	// This method flip the terminals of ASP GLU PHY ARG and TYR.  
	public static void flipTerm(Residue res) {
		int restype = res.type;
		if ((restype==ASP) || (restype==GLU) || (restype==PHE) || (restype==ARG) ||  (restype==TYR)) {
			Atom a1=null,a2=null;
			double tx,ty,tz;
			if (restype==ASP) {//ASP
				a1 = res.atoms().getAtom("OD1");
				a2 = res.atoms().getAtom("OD2");
			}   
			if (restype==GLU) {//GLU
				a1 = res.atoms().getAtom("OE1");
				a2 = res.atoms().getAtom("OE2");
			}   
			if (restype==PHE) {//PHE
				a1 = res.atoms().getAtom("CD1");
				a2 = res.atoms().getAtom("CD2");
			}   
			if (restype==ARG) {//ARG
				a1 = res.atoms().getAtom("NH1");
				a2 = res.atoms().getAtom("NH2");
			}   
			if (restype==TYR) {//TYR
				a1 = res.atoms().getAtom("CD1");
				a2 = res.atoms().getAtom("CD2");
			}   
			tx = a1.x();
			ty = a1.y();
			tz = a1.z();
			a1.setXYZ(a2.x(),a2.y(),a2.z());
			a2.setXYZ(tx,ty,tz);
			if (restype==PHE) {//PHE
				a1 = res.atoms().getAtom("CE1");
				a2 = res.atoms().getAtom("CE2");
			}   
			if (restype==TYR) {//TYR
				a1 = res.atoms().getAtom("CE1");
				a2 = res.atoms().getAtom("CE2");
			}   
			tx = a1.x();
			ty = a1.y();
			tz = a1.z();
			a1.setXYZ(a2.x(),a2.y(),a2.z());
			a2.setXYZ(tx,ty,tz);
		}
	} 


	/**
	 * This method correct all the sidechains in the protein with flat termini (i.e. ASP,GLU,PHE etc.) whose abs(chi_terminal) value is above 
	 * 90 degrees. The correction is by flipping the terminal. This conforms with the IUPAC convention.
	 */
	public static void correctNomenclature(Protein prot) {
		double[][] chis = getSidechainTorsions(prot);
		for (int r=0; r<prot.residues().size() ; r++){
			int resType=prot.residues().residueAt(r).type;
			int resNum=prot.residues().residueAt(r).number;
			if ((resType==ASP) && // ASP
					(Math.abs(chis[resNum][1])>(Math.PI/2)))
				flipTerm(prot.residues().residueAt(r));
			if ((resType==GLU) && // GLU
					(Math.abs(chis[resNum][2])>(Math.PI/2)))
				flipTerm(prot.residues().residueAt(r));
			if ((resType==PHE) && // PHE
					(Math.abs(chis[resNum][1])>(Math.PI/2)))
				flipTerm(prot.residues().residueAt(r));
			if ((resType==TYR) && // TYR
					(Math.abs(chis[resNum][1])>(Math.PI/2)))
				flipTerm(prot.residues().residueAt(r));
			if (resType==ARG) {  // ARG
				Atom cd = prot.residues().residueAt(r).atoms().getAtom("CD");
				Atom nh1 = prot.residues().residueAt(r).atoms().getAtom("NH1");
				Atom nh2 = prot.residues().residueAt(r).atoms().getAtom("NH2");
				if ((cd!=null)&&(nh1!=null)&&(nh2!=null)&&(cd.distanceFrom(nh2) < cd.distanceFrom(nh1)))
					flipTerm(prot.residues().residueAt(r));				
			}
		}
	}
	
	
	// This method puts all the residues of prot, that have at least ONE NON-FROZEN atom, into a random rotamer from the library.
	// It also returns an Nx4 double array where row n is the phi,psi,type,{rotamer prob} of residue n.
	// Note that the distance matrix is discarded after this method is run, as it is the best way to ensure,
	// that it is updated properly after this sidechain change.
	public static double[][] putIntoRandomRot(Protein prot, DistanceMatrix dm, DunbrackLib lib,double[] randomVals) { 

		double[][] pp = phipsi(prot,dm);
		double[] tmprot;

		int maxRes = -1;
		for (int kk=0 ; kk<prot.atoms().size() ; kk++) 
			if (prot.atoms().atomAt(kk).residueNumber() > maxRes)
				maxRes = prot.atoms().atomAt(kk).residueNumber();
		double[][] largepp = new double[maxRes+1][];
		for (int kk=0 ; kk<largepp.length ; kk++) 
			if ((prot.residue(kk)!=null) && (prot.residue(kk).type>-1))
				if ((kk<pp.length) && (pp[kk]!=null))
					largepp[kk] = pp[kk];
				else {
					largepp[kk] = new double[3];
					largepp[kk][0] = -60*Math.PI/180.0;
					largepp[kk][1] = -40*Math.PI/180.0;
					largepp[kk][2] = prot.residue(kk).type;
				}
		pp = largepp;		

		// putting Rot1 into all instances of this residue type
		double[][] tmppp = new double[pp.length][];
		for (int kk=0 ; kk<pp.length ; kk++) {
			if ((pp[kk]!=null) && (pp[kk][2]>-1) && (prot.residue(kk).atoms().defrostedAtoms().size()>0)) {
				if (pp[kk][2] ==ALA)
					ResidueBuilder.build(prot.residue(kk),0,null);
				else if (pp[kk][2]!=GLY) {
					int rotNum = (int) (lib.getRotamerNum((int) pp[kk][2] , pp[kk][0] , pp[kk][1])
							*randomVals[kk%randomVals.length]);
					tmprot = lib.getRotamer((int) pp[kk][2] , pp[kk][0] , pp[kk][1] , rotNum);
					ResidueBuilder.build(prot.residue(kk),prot.residue(kk).type,tmprot);
					tmppp[kk] = new double[4];
					tmppp[kk][0] = pp[kk][0];
					tmppp[kk][1] = pp[kk][1];
					tmppp[kk][2] = pp[kk][2];
					tmppp[kk][3] = lib.getRotamerProb((int) pp[kk][2] , pp[kk][0] , pp[kk][1] , rotNum); 	
				}
			}
		}

		// Discarding the distance matrix.
		dm.terminator.kill("The distance matrix is not updated correctly in PutIntoRot1. It has to be created anew.");
		return tmppp;    
	}  // of putIntoRandomRot

	
	/**
	 *  Return the coordinates [x,y,z] of the representative atom of residue 'res'.
	 *  
	 *  If the coordinates cannot be found (like in ALA) a runtime error is thrown.
	 *  
	 **/
	public static double[] findRepCoors(Residue res) {
		double[] coors = new double[3];
		Atom atom1=null;
		Atom atom2 = null;
		switch (res.type) {
		case 1:  // CYS
			atom1 = res.atoms().findAtomInList("SG", res.number);
			break;
		case 2:  // ASP
			atom1 = res.atoms().findAtomInList("OD1", res.number);
			atom2 = res.atoms().findAtomInList("OD2", res.number);
			break;
		case 3:  // GLU
			atom1 = res.atoms().findAtomInList("OE1", res.number);
			atom2 = res.atoms().findAtomInList("OE2", res.number);
			break;
		case 4:  // PHE
			atom1 = res.atoms().findAtomInList("CZ", res.number);
			break;
		case 6:  // HIS
			atom1 = res.atoms().findAtomInList("CE1", res.number);
			atom2 = res.atoms().findAtomInList("NE2", res.number);
			break;
		case 7:  // ILE
			atom1 = res.atoms().findAtomInList("CD1", res.number);
			break;
		case 8:  // LYS
			atom1 = res.atoms().findAtomInList("NZ", res.number);
			break;
		case 9:  // LEU
			atom1 = res.atoms().findAtomInList("CD1", res.number);
			atom2 = res.atoms().findAtomInList("CD2", res.number);
			break;
		case 10:  // MET
			atom1 = res.atoms().findAtomInList("CE", res.number);
			break;
		case 11:  // ASN
			atom1 = res.atoms().findAtomInList("OD1", res.number);
			atom2 = res.atoms().findAtomInList("ND2", res.number);
			break;
		case 12:  // PRO
			atom1 = res.atoms().findAtomInList("CG", res.number);
			break;
		case 13:  // GLN
			atom1 = res.atoms().findAtomInList("OE1", res.number);
			atom2 = res.atoms().findAtomInList("NE2", res.number);
			break;
		case 14:  // ARG
			atom1 = res.atoms().findAtomInList("NH1", res.number);
			atom2 = res.atoms().findAtomInList("NH2", res.number);
			break;
		case 15:  // SER
			atom1 = res.atoms().findAtomInList("OG", res.number);
			break;
		case 16:  // THR
			atom1 = res.atoms().findAtomInList("OG1", res.number);
			break;
		case 17:  // VAL
			atom1 = res.atoms().findAtomInList("CG1", res.number);
			atom2 = res.atoms().findAtomInList("CG2", res.number);
			break;
		case 18:  // TRP
			atom1 = res.atoms().findAtomInList("CH2", res.number);
			break;
		case 19:  // TYR
			atom1 = res.atoms().findAtomInList("OH", res.number);
			break;
		default:
			throw new RuntimeException("A residue type with no sidechain");				
		}
		if (atom1==null)
			throw new RuntimeException("Could not find representative atom.");				
		if (atom2==null) {
			coors[0] = atom1.x();
			coors[1] = atom1.y();
			coors[2] = atom1.z();
		} else {
			coors[0] = 0.5*(atom1.x() + atom2.x());
			coors[1] = 0.5*(atom1.y() + atom2.y());
			coors[2] = 0.5*(atom1.z() + atom2.z());			
		}
		return coors;
	}
	
	/**
	 * Return true is the residue 'res' has all the atoms in its sidechain.
	 */
	public static boolean isSidechainComplete(Residue res) {
		if (res==null)
			return true;
		AtomList tmpList = res.atoms().sidechain();
		switch (res.type) {
		case 0:  // ALA
			return true;
		case 1:  // CYS
			if (tmpList.size()==1)
				return true;
			else
				return false;
		case 2:  // ASP
			if (tmpList.size()==3)
				return true;
			else
				return false;
		case 3:  // GLU
			if (tmpList.size()==4)
				return true;
			else
				return false;
		case 4:  // PHE
			if (tmpList.size()==6)
				return true;
			else
				return false;
		case 6:  // HIS
			if (tmpList.size()==7)
				return true;
			else
				return false;
		case 7:  // ILE
			if (tmpList.size()==3)
				return true;
			else
				return false;
		case 8:  // LYS
			if (tmpList.size()==4)
				return true;
			else
				return false;
		case 9:  // LEU
			if (tmpList.size()==3)
				return true;
			else
				return false;
		case 10:  // MET
			if (tmpList.size()==3)
				return true;
			else
				return false;
		case 11:  // ASN
			if (tmpList.size()==5)
				return true;
			else
				return false;
		case 12:  // PRO
			if (tmpList.size()==2)
				return true;
			else
				return false;
		case 13:  // GLN
			if (tmpList.size()==6)
				return true;
			else
				return false;
		case 14:  // ARG
			if (tmpList.size()==7)
				return true;
			else
				return false;
		case 15:  // SER
			if (tmpList.size()==1)
				return true;
			else
				return false;
		case 16:  // THR
			if (tmpList.size()==2)
				return true;
			else
				return false;
		case 17:  // VAL
			if (tmpList.size()==2)
				return true;
			else
				return false;
		case 18:  // TRP
			if (tmpList.size()==10)
				return true;
			else
				return false;
		case 19:  // TYR
			if (tmpList.size()==7)
				return true;
			else
				return false;
		default:
			return true;
		}		
	}
	
	/**
	 *  This method assumes that the protein is in ROT1 representaion.
	 *  
	 *  It puts the CB in a representative locus for each sidechain. 
	 *   
	 **/
	public static void putCBinSCcenter(Residue res) {
		Atom atom1=null;
		Atom atom2 = null;
		switch (res.type) {
		case 0:  // ALA
			return;
		case 1:  // CYS
			atom1 = res.atoms().findAtomInList("SG", res.number);
			break;
		case 2:  // ASP
			atom1 = res.atoms().findAtomInList("OD1", res.number);
			atom2 = res.atoms().findAtomInList("OD2", res.number);
			break;
		case 3:  // GLU
			atom1 = res.atoms().findAtomInList("OE1", res.number);
			atom2 = res.atoms().findAtomInList("OE2", res.number);
			break;
		case 4:  // PHE
			atom1 = res.atoms().findAtomInList("CZ", res.number);
			atom2 = res.atoms().findAtomInList("CG", res.number);
			break;
		case 5:  // GLY
			return;
		case 6:  // HIS
			atom1 = res.atoms().findAtomInList("CE1", res.number);
			atom2 = res.atoms().findAtomInList("NE2", res.number);
			break;
		case 7:  // ILE
			atom1 = res.atoms().findAtomInList("CD1", res.number);
			break;
		case 8:  // LYS
			atom1 = res.atoms().findAtomInList("NZ", res.number);
			break;
		case 9:  // LEU
			atom1 = res.atoms().findAtomInList("CD1", res.number);
			atom2 = res.atoms().findAtomInList("CD2", res.number);
			break;
		case 10:  // MET
			atom1 = res.atoms().findAtomInList("CE", res.number);
			break;
		case 11:  // ASN
			atom1 = res.atoms().findAtomInList("OD1", res.number);
			atom2 = res.atoms().findAtomInList("ND2", res.number);
			break;
		case 12:  // PRO
			atom1 = res.atoms().findAtomInList("CG", res.number);
			break;
		case 13:  // GLN
			atom1 = res.atoms().findAtomInList("OE1", res.number);
			atom2 = res.atoms().findAtomInList("NE2", res.number);
			break;
		case 14:  // ARG
			atom1 = res.atoms().findAtomInList("CZ", res.number);
			break;
		case 15:  // SER
			atom1 = res.atoms().findAtomInList("OG", res.number);
			break;
		case 16:  // THR
			atom1 = res.atoms().findAtomInList("OG1", res.number);
			atom2 = res.atoms().findAtomInList("CG2", res.number);
			break;
		case 17:  // VAL
			atom1 = res.atoms().findAtomInList("CG1", res.number);
			atom2 = res.atoms().findAtomInList("CG2", res.number);
			break;
		case 18:  // TRP
			atom1 = res.atoms().findAtomInList("CH2", res.number);
			atom2 = res.atoms().findAtomInList("CD2", res.number);
			break;
		case 19:  // TYR
			atom1 = res.atoms().findAtomInList("OH", res.number);
			atom2 = res.atoms().findAtomInList("CG", res.number);
			break;
		default:
			throw new RuntimeException("A residue type with no sidechain");				
		}
		if (atom1==null)
			throw new RuntimeException("Could not find representative atom.");				
		if (atom2==null) {
			res.atoms().findAtomInList("CB", res.number).setXYZ(atom1.x()+0.05, atom1.y()+0.05, atom1.z()+0.05);
		} else {
			res.atoms().findAtomInList("CB", res.number).setXYZ(0.5*(atom1.x()+atom2.x())+0.05, 0.5*(atom1.y()+atom2.y())+0.05, 0.5*(atom1.z()+atom2.z())+0.05);
		}
	}

	public static void putOinSCcenter(Residue res) {
		Atom atom1=null;
		Atom atom2 = null;
		switch (res.type) {
		case 0:  // ALA
			atom1 = res.atoms().findAtomInList("CB", res.number);
			break;
		case 1:  // CYS
			atom1 = res.atoms().findAtomInList("SG", res.number);
			break;
		case 2:  // ASP
			atom1 = res.atoms().findAtomInList("OD1", res.number);
			atom2 = res.atoms().findAtomInList("OD2", res.number);
			break;
		case 3:  // GLU
			atom1 = res.atoms().findAtomInList("OE1", res.number);
			atom2 = res.atoms().findAtomInList("OE2", res.number);
			break;
		case 4:  // PHE
			atom1 = res.atoms().findAtomInList("CZ", res.number);
			atom2 = res.atoms().findAtomInList("CG", res.number);
			break;
		case 5:  // GLY
			return;
		case 6:  // HIS
			atom1 = res.atoms().findAtomInList("CE1", res.number);
			atom2 = res.atoms().findAtomInList("NE2", res.number);
			break;
		case 7:  // ILE
			atom1 = res.atoms().findAtomInList("CD1", res.number);
			break;
		case 8:  // LYS
			atom1 = res.atoms().findAtomInList("NZ", res.number);
			break;
		case 9:  // LEU
			atom1 = res.atoms().findAtomInList("CD1", res.number);
			atom2 = res.atoms().findAtomInList("CD2", res.number);
			break;
		case 10:  // MET
			atom1 = res.atoms().findAtomInList("CE", res.number);
			break;
		case 11:  // ASN
			atom1 = res.atoms().findAtomInList("OD1", res.number);
			atom2 = res.atoms().findAtomInList("ND2", res.number);
			break;
		case 12:  // PRO
			atom1 = res.atoms().findAtomInList("CG", res.number);
			break;
		case 13:  // GLN
			atom1 = res.atoms().findAtomInList("OE1", res.number);
			atom2 = res.atoms().findAtomInList("NE2", res.number);
			break;
		case 14:  // ARG
			atom1 = res.atoms().findAtomInList("CZ", res.number);
			break;
		case 15:  // SER
			atom1 = res.atoms().findAtomInList("OG", res.number);
			break;
		case 16:  // THR
			atom1 = res.atoms().findAtomInList("OG1", res.number);
			atom2 = res.atoms().findAtomInList("CG2", res.number);
			break;
		case 17:  // VAL
			atom1 = res.atoms().findAtomInList("CG1", res.number);
			atom2 = res.atoms().findAtomInList("CG2", res.number);
			break;
		case 18:  // TRP
			atom1 = res.atoms().findAtomInList("CH2", res.number);
			atom2 = res.atoms().findAtomInList("CD2", res.number);
			break;
		case 19:  // TYR
			atom1 = res.atoms().findAtomInList("OH", res.number);
			atom2 = res.atoms().findAtomInList("CG", res.number);
			break;
		default:
			throw new RuntimeException("A residue type with no sidechain");				
		}
		if (atom1==null)
			throw new RuntimeException("Could not find representative atom.");				
		if (atom2==null) {
			res.atoms().findAtomInList("O", res.number).setXYZ(atom1.x()+0.05, atom1.y()+0.05, atom1.z()+0.05);
		} else {
			res.atoms().findAtomInList("O", res.number).setXYZ(0.5*(atom1.x()+atom2.x())+0.05, 0.5*(atom1.y()+atom2.y())+0.05, 0.5*(atom1.z()+atom2.z())+0.05);
		}
	}	

} // of class

