package meshi.geometry;

import meshi.molecularElements.AtomList;
import meshi.molecularElements.Protein;
import meshi.util.mathTools.Jama.Matrix;


public class Helixer {

AtomList full = new AtomList();
AtomList bb = new AtomList();
double cmx=0,cmy=0,cmz=0;
double[][] rotBack;
int nbb,nfull;
double[][] bbcoor;
double[][] fullcoor;


public Helixer(Protein prot , int resBegin , int resEnd) {

	for (int c=0 ; c<prot.atoms().size() ; c++) 
		if ((prot.atoms().atomAt(c).residueNumber()>=resBegin) &&
			(prot.atoms().atomAt(c).residueNumber()<=resEnd)) {
			full.add(prot.atoms().atomAt(c));
			if (prot.atoms().atomAt(c).isBackbone) { 
				bb.add(prot.atoms().atomAt(c));
				cmx += prot.atoms().atomAt(c).x();
				cmy += prot.atoms().atomAt(c).y();
				cmz += prot.atoms().atomAt(c).z();
			}
	}
	nbb = bb.size();
	nfull = full.size();
	cmx /= nbb;
	cmy /= nbb;
	cmz /= nbb;
	
	bbcoor = new double[nbb][3];
	fullcoor = new double[nfull][3];
	for (int c=0 ; c<nbb ; c++) {
		bbcoor[c][0] = bb.atomAt(c).x()-cmx;
		bbcoor[c][1] = bb.atomAt(c).y()-cmy;
		bbcoor[c][2] = bb.atomAt(c).z()-cmz;
	} 
	for (int c=0 ; c<nfull ; c++) {
		fullcoor[c][0] = full.atomAt(c).x()-cmx;
		fullcoor[c][1] = full.atomAt(c).y()-cmy;
		fullcoor[c][2] = full.atomAt(c).z()-cmz;
	} 
	
	double[][] m  = new double[3][3]; // Covarient matrix
	for (int c=0 ; c<nbb ; c++) {
		m[0][0] += bbcoor[c][0]*bbcoor[c][0];
		m[0][1] += bbcoor[c][0]*bbcoor[c][1];
		m[0][2] += bbcoor[c][0]*bbcoor[c][2];
		m[1][0] = m[0][1];
		m[1][1] += bbcoor[c][1]*bbcoor[c][1];
		m[1][2] += bbcoor[c][1]*bbcoor[c][2];
		m[2][0] = m[0][2];
		m[2][1] = m[1][2];
		m[2][2] += bbcoor[c][2]*bbcoor[c][2];
	}
	
	Matrix mm = new Matrix(m);
	double[] lambda = mm.eig().getRealEigenvalues();
	Matrix eigv = mm.eig().getV();
	int maxind = -1;
	if ((Math.abs(lambda[0])>=Math.abs(lambda[1])) && (Math.abs(lambda[0])>=Math.abs(lambda[2])))
		maxind = 0;
	else if ((Math.abs(lambda[1])>=Math.abs(lambda[0])) && (Math.abs(lambda[1])>=Math.abs(lambda[2])))
			maxind = 1;
		else
			maxind = 2;
	
	double[] pca1 = new double[3];

	pca1[0] = eigv.get(0,maxind);
	pca1[1] = eigv.get(1,maxind);
	pca1[2] = eigv.get(2,maxind);
	double pca1_norm = norm(pca1);
	pca1[0] /= pca1_norm;
	pca1[1] /= pca1_norm;
	pca1[2] /= pca1_norm;
	
	// fliping the pca, so that it is closest to the helix end
	if (dis(pca1[0],pca1[1],pca1[2],
		prot.residue(resEnd).ca().x(),
		prot.residue(resEnd).ca().y(),
		prot.residue(resEnd).ca().z()) >  dis(-pca1[0],-pca1[1],-pca1[2],
		prot.residue(resEnd).ca().x(),
		prot.residue(resEnd).ca().y(),
		prot.residue(resEnd).ca().z())) {
			pca1[0] = -pca1[0];
			pca1[1] = -pca1[1];
			pca1[2] = -pca1[2];
	}

	// The rotation axis
	double rx = pca1[1];
	double ry = -pca1[0];
	rx /= norm(pca1[0],pca1[1],0);
	ry /= norm(pca1[0],pca1[1],0);
	
	// The angle to the z
	double ang = Math.acos(pca1[2]);
	double cosang = pca1[2];
	double sinang = Math.sin(ang);
	
	
	double[][] rot = new double[3][3];
	rot[0][0] = cosang + (1-cosang)*rx*rx;
	rot[0][1] = (1-cosang)*rx*ry;
	rot[0][2] = sinang*ry;
	rot[1][0] = (1-cosang)*rx*ry;	
	rot[1][1] = cosang + (1-cosang)*ry*ry;
	rot[1][2] = -sinang*rx;
	rot[2][0] = -sinang*ry;
	rot[2][1] = sinang*rx;
	rot[2][2] = cosang;
	
	// Rotating:
	for (int c=0 ; c<nbb ; c++) {
		bbcoor[c] = matTimesVec(rot,bbcoor[c]);
	} 
	for (int c=0 ; c<nfull ; c++) {
		fullcoor[c] = matTimesVec(rot,fullcoor[c]);
	}
	
	
	// The inverse rotation matrix we keep:
	ang = -ang;
	cosang = pca1[2];
	sinang = Math.sin(ang);
	
	
	rotBack = new double[3][3];
	rotBack[0][0] = cosang + (1-cosang)*rx*rx;
	rotBack[0][1] = (1-cosang)*rx*ry;
	rotBack[0][2] = sinang*ry;
	rotBack[1][0] = (1-cosang)*rx*ry;	
	rotBack[1][1] = cosang + (1-cosang)*ry*ry;
	rotBack[1][2] = -sinang*rx;
	rotBack[2][0] = -sinang*ry;
	rotBack[2][1] = sinang*rx;
	rotBack[2][2] = cosang;
	 

/*  System.out.println();
	for (int c=0 ; c<nbb ; c++) {
		System.out.println(bbcoor[c][0] + " " + bbcoor[c][1] + " " + bbcoor[c][2]);
		System.out.println(c + " " + fullcoor[c][0] + " " + fullcoor[c][1] + " " + fullcoor[c][2]);
		System.out.println(c + " " + tmp[0] + " " + tmp[1] + " " + tmp[2]);
	}
*/ 
}	// Constructor
	
public void rockNroll(double angx, double angy, double angz, double offsetx, double offsety, double offsetz) {
	double[][] zrotmat = {{Math.cos(angz),Math.sin(angz),0},
						{-Math.sin(angz),Math.cos(angz),0},
						{0,0,1}};	
	double[][] yrotmat = {{Math.cos(angy),0,Math.sin(angy)},
						{0,1,0},	
						{-Math.sin(angy),0,Math.cos(angy)}};	
	double[][] xrotmat = {{1,0,0},
						{0,Math.cos(angx),Math.sin(angx)},
						{0,-Math.sin(angx),Math.cos(angx)}};	

	for (int c=0 ; c<nfull ; c++) {
		double[] tmp = matTimesVec(rotBack,matTimesVec(xrotmat,matTimesVec(yrotmat,matTimesVec(zrotmat,fullcoor[c]))));
		full.atomAt(c).setX(tmp[0]+cmx+offsetx*rotBack[0][0]+offsety*rotBack[0][1]+offsetz*rotBack[0][2]);
		full.atomAt(c).setY(tmp[1]+cmy+offsetx*rotBack[1][0]+offsety*rotBack[1][1]+offsetz*rotBack[1][2]);
		full.atomAt(c).setZ(tmp[2]+cmz+offsetx*rotBack[2][0]+offsety*rotBack[2][1]+offsetz*rotBack[2][2]);
	} 	
}

	
	

public static double dis(double x1,double y1,double z1,double x2,double y2,double z2) {
	return Math.sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}

	
public static double norm(double x,double y,double z) {
	double[] tmp = {x,y,z};
	return norm(tmp);
}


public static double norm(double[] ar) {
	double sum = 0.0;
	for (int c=0 ; c<ar.length ; c++) {
		sum += ar[c]*ar[c];
	}
	return Math.sqrt(sum);
}	

public static double[] matTimesVec(double[][] m, double[] v) {
	double[] result = new double[v.length];
	for (int c=0 ; c<v.length ; c++)
		for (int d=0 ; d<v.length ; d++)
			result[c] += m[c][d]*v[d];
	return result; 
}
	
}
