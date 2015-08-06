package meshi.applications.TriC;

import meshi.applications.corpus.Corpus;

/** 
 * This class calculates the sequence similarity between subunits.
 * 
 * @author Nir
 *
 */

public class Calculate_Sequence_Similarity {
	
	
	public static void main(String[] args) {
		String conversion = "ABGDEZHQ";
		String conversion2 = "ACDEFGHIKLMNPQRSTVWY";
		TricAlignment tricAlignment = new TricAlignment();
		VectorOfOrganisms vec = new VectorOfOrganisms();
		vec.add(new ReadTriCofOrganism());
		vec.add(new ReadTriCofOrganism("human.txt"));
		vec.add(new ReadTriCofOrganism("zebrafish.txt"));
		vec.add(new ReadTriCofOrganism("ciona.txt"));
		vec.add(new ReadTriCofOrganism("drosoph.txt"));
		vec.add(new ReadTriCofOrganism("Celegance.txt"));
		vec.add(new ReadTriCofOrganism("arabidopsis.txt"));
		vec.add(new ReadTriCofOrganism("yeast.txt"));
		vec.add(new ReadTriCofOrganism("neospora.txt"));
		vec.add(new ReadTriCofOrganism("candida.txt"));
		vec.add(new ReadTriCofOrganism("mold.txt"));
		vec.add(new ReadTriCofOrganism("plasmodium.txt"));
		vec.add(new ReadTriCofOrganism("paramecium.txt"));
		
		double[] count = new double[8];
		double[] sum = new double[8];
		int alignmentLength = 555;
		for (int pos=0 ; pos<alignmentLength ; pos++) {
			// Check if this is non-gap
			boolean foundGap = false;
			for (int unit=0 ; unit<8 ; unit++) {
				String tmp = vec.getProfile(conversion.charAt(unit), pos);
				if (tmp.indexOf("-")!=-1) {
					foundGap = true;
				}				
			}
			if (!foundGap) {
				for (int unit1=0 ; unit1<8 ; unit1++) {
					for (int unit2=0 ; unit2<8 ; unit2++) {
						if (unit1!=unit2) {
							String prof1 = vec.getProfile(conversion.charAt(unit1), pos);
							String prof2 = vec.getProfile(conversion.charAt(unit2), pos);
							for (int c1=0 ; c1<prof1.length() ; c1++) {
								for (int c2=0 ; c2<prof1.length() ; c2++) {
									count[unit1]++;
//									sum[unit1] += Corpus.blosum62[conversion2.indexOf(""+prof1.charAt(c1))][conversion2.indexOf(""+prof2.charAt(c2))]; // Sum of blossum scores.

									if (Corpus.blosum62[conversion2.indexOf(""+prof1.charAt(c1))][conversion2.indexOf(""+prof2.charAt(c2))] > 0) {
										sum[unit1]++;
									}								
									
									
//									if (prof1.charAt(c1)==prof2.charAt(c2)) {
//										sum[unit1]++;
//									}								
								}
							}							
						}						
					}					
				}
			}
		}
		
		for (int unit1=0 ; unit1<8 ; unit1++) {
			System.out.println(conversion.charAt(unit1)+": "+(sum[unit1]/count[unit1]) + "    " + count[unit1]);
		}
	}

}
