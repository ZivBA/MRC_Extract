package meshi.applications.TriC;

public class GetTwoPositions {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		ReadTriCofOrganism yeast = new ReadTriCofOrganism("yeast.txt");
		yeast.printAlign();
		System.exit(0);
		
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
		vec.add(new ReadTriCofOrganism("trypansoma.txt"));
		vec.add(new ReadTriCofOrganism("leishmania.txt"));
		vec.add(new ReadTriCofOrganism("tricto.txt"));
		
//        // Interaction center 1
//		System.out.println(vec.getProfile('A', 41));
//		System.out.println(vec.getProfile('B', 41));
//		System.out.println(vec.getProfile('G', 41));
//		System.out.println(vec.getProfile('D', 41));
//		System.out.println(vec.getProfile('E', 41));
//		System.out.println(vec.getProfile('H', 41));
//		System.out.println(vec.getProfile('Q', 41));
//		System.out.println(vec.getProfile('Z', 41));
//		System.out.println();
//		System.out.println(vec.getProfile('A', 525));
//		System.out.println(vec.getProfile('B', 525));
//		System.out.println(vec.getProfile('G', 525));
//		System.out.println(vec.getProfile('D', 525));
//		System.out.println(vec.getProfile('E', 525));
//		System.out.println(vec.getProfile('H', 525));
//		System.out.println(vec.getProfile('Q', 525));
//		System.out.println(vec.getProfile('Z', 525));
		
//        // Interaction center 2
//		System.out.println(vec.getProfile('A', 51));
//		System.out.println(vec.getProfile('B', 51));
//		System.out.println(vec.getProfile('G', 51));
//		System.out.println(vec.getProfile('D', 51));
//		System.out.println(vec.getProfile('E', 51));
//		System.out.println(vec.getProfile('H', 51));
//		System.out.println(vec.getProfile('Q', 51));
//		System.out.println(vec.getProfile('Z', 51));
//		System.out.println();
//		System.out.println(vec.getProfile('A', 50));
//		System.out.println(vec.getProfile('B', 50));
//		System.out.println(vec.getProfile('G', 50));
//		System.out.println(vec.getProfile('D', 50));
//		System.out.println(vec.getProfile('E', 50));
//		System.out.println(vec.getProfile('H', 50));
//		System.out.println(vec.getProfile('Q', 50));
//		System.out.println(vec.getProfile('Z', 50));
//		System.out.println();
//		System.out.println(vec.getProfile('A', 72));
//		System.out.println(vec.getProfile('B', 72));
//		System.out.println(vec.getProfile('G', 72));
//		System.out.println(vec.getProfile('D', 72));
//		System.out.println(vec.getProfile('E', 72));
//		System.out.println(vec.getProfile('H', 72));
//		System.out.println(vec.getProfile('Q', 72));
//		System.out.println(vec.getProfile('Z', 72));
//		System.out.println();
//		System.out.println(vec.getProfile('A', 69));
//		System.out.println(vec.getProfile('B', 69));
//		System.out.println(vec.getProfile('G', 69));
//		System.out.println(vec.getProfile('D', 69));
//		System.out.println(vec.getProfile('E', 69));
//		System.out.println(vec.getProfile('H', 69));
//		System.out.println(vec.getProfile('Q', 69));
//		System.out.println(vec.getProfile('Z', 69));
		
//		// The mess around the histedine - We will get back to it latter
//	    int res = 380;
//		System.out.println(vec.getProfile('A', res));
//		System.out.println(vec.getProfile('B', res));
//		System.out.println(vec.getProfile('G', res));
//		System.out.println(vec.getProfile('D', res));
//		System.out.println(vec.getProfile('E', res));
//		System.out.println(vec.getProfile('H', res));
//		System.out.println(vec.getProfile('Q', res));
//		System.out.println(vec.getProfile('Z', res));
//		System.out.println();
//
//		res = 382;
//		System.out.println(vec.getProfile('A', res));
//		System.out.println(vec.getProfile('B', res));
//		System.out.println(vec.getProfile('G', res));
//		System.out.println(vec.getProfile('D', res));
//		System.out.println(vec.getProfile('E', res));
//		System.out.println(vec.getProfile('H', res));
//		System.out.println(vec.getProfile('Q', res));
//		System.out.println(vec.getProfile('Z', res));
//		System.out.println();
//		
//		res = 73;
//		System.out.println(vec.getProfile('A', res));
//		System.out.println(vec.getProfile('B', res));
//		System.out.println(vec.getProfile('G', res));
//		System.out.println(vec.getProfile('D', res));
//		System.out.println(vec.getProfile('E', res));
//		System.out.println(vec.getProfile('H', res));
//		System.out.println(vec.getProfile('Q', res));
//		System.out.println(vec.getProfile('Z', res));
//		System.out.println();
//		
//		res = 76;
//		System.out.println(vec.getProfile('A', res));
//		System.out.println(vec.getProfile('B', res));
//		System.out.println(vec.getProfile('G', res));
//		System.out.println(vec.getProfile('D', res));
//		System.out.println(vec.getProfile('E', res));
//		System.out.println(vec.getProfile('H', res));
//		System.out.println(vec.getProfile('Q', res));
//		System.out.println(vec.getProfile('Z', res));
//		System.out.println();

			
		int resInThermoB = 333;
		System.out.println("Therm " + resInThermoB + " " + tricAlignment.getAAinSeq('J', resInThermoB) +
				" " + tricAlignment.getAAinSeq('J', resInThermoB,'I'));
		System.out.println("unitA " + tricAlignment.getAAinSeq('J', resInThermoB,'A') + 
				vec.getProfile('A', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitB " + tricAlignment.getAAinSeq('J', resInThermoB,'B') + 
				vec.getProfile('B', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitG " + tricAlignment.getAAinSeq('J', resInThermoB,'G') + 
				vec.getProfile('G', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitD " + tricAlignment.getAAinSeq('J', resInThermoB,'D') + 
				vec.getProfile('D', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitE " + tricAlignment.getAAinSeq('J', resInThermoB,'E') + 
				vec.getProfile('E', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitH " + tricAlignment.getAAinSeq('J', resInThermoB,'H') + 
				vec.getProfile('H', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitQ " + tricAlignment.getAAinSeq('J', resInThermoB,'Q') + 
				vec.getProfile('Q', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitZ " + tricAlignment.getAAinSeq('J', resInThermoB,'Z') + 
				vec.getProfile('Z', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println();		

		resInThermoB = 210;
		System.out.println("Therm " + resInThermoB + " " + tricAlignment.getAAinSeq('J', resInThermoB) +
				" " + tricAlignment.getAAinSeq('J', resInThermoB,'I'));
		System.out.println("unitA " + tricAlignment.getAAinSeq('J', resInThermoB,'A') + 
				vec.getProfile('A', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitB " + tricAlignment.getAAinSeq('J', resInThermoB,'B') + 
				vec.getProfile('B', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitG " + tricAlignment.getAAinSeq('J', resInThermoB,'G') + 
				vec.getProfile('G', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitD " + tricAlignment.getAAinSeq('J', resInThermoB,'D') + 
				vec.getProfile('D', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitE " + tricAlignment.getAAinSeq('J', resInThermoB,'E') + 
				vec.getProfile('E', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitH " + tricAlignment.getAAinSeq('J', resInThermoB,'H') + 
				vec.getProfile('H', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitQ " + tricAlignment.getAAinSeq('J', resInThermoB,'Q') + 
				vec.getProfile('Q', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitZ " + tricAlignment.getAAinSeq('J', resInThermoB,'Z') + 
				vec.getProfile('Z', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println();		

		resInThermoB = 209;
		System.out.println("Therm " + resInThermoB + " " + tricAlignment.getAAinSeq('J', resInThermoB) +
				" " + tricAlignment.getAAinSeq('J', resInThermoB,'I'));
		System.out.println("unitA " + tricAlignment.getAAinSeq('J', resInThermoB,'A') + 
				vec.getProfile('A', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitB " + tricAlignment.getAAinSeq('J', resInThermoB,'B') + 
				vec.getProfile('B', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitG " + tricAlignment.getAAinSeq('J', resInThermoB,'G') + 
				vec.getProfile('G', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitD " + tricAlignment.getAAinSeq('J', resInThermoB,'D') + 
				vec.getProfile('D', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitE " + tricAlignment.getAAinSeq('J', resInThermoB,'E') + 
				vec.getProfile('E', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitH " + tricAlignment.getAAinSeq('J', resInThermoB,'H') + 
				vec.getProfile('H', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitQ " + tricAlignment.getAAinSeq('J', resInThermoB,'Q') + 
				vec.getProfile('Q', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitZ " + tricAlignment.getAAinSeq('J', resInThermoB,'Z') + 
				vec.getProfile('Z', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println();	

		resInThermoB = 208;
		System.out.println("Therm " + resInThermoB + " " + tricAlignment.getAAinSeq('J', resInThermoB) +
				" " + tricAlignment.getAAinSeq('J', resInThermoB,'I'));
		System.out.println("unitA " + tricAlignment.getAAinSeq('J', resInThermoB,'A') + 
				vec.getProfile('A', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitB " + tricAlignment.getAAinSeq('J', resInThermoB,'B') + 
				vec.getProfile('B', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitG " + tricAlignment.getAAinSeq('J', resInThermoB,'G') + 
				vec.getProfile('G', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitD " + tricAlignment.getAAinSeq('J', resInThermoB,'D') + 
				vec.getProfile('D', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitE " + tricAlignment.getAAinSeq('J', resInThermoB,'E') + 
				vec.getProfile('E', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitH " + tricAlignment.getAAinSeq('J', resInThermoB,'H') + 
				vec.getProfile('H', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitQ " + tricAlignment.getAAinSeq('J', resInThermoB,'Q') + 
				vec.getProfile('Q', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitZ " + tricAlignment.getAAinSeq('J', resInThermoB,'Z') + 
				vec.getProfile('Z', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println();	

		resInThermoB = 207;
		System.out.println("Therm " + resInThermoB + " " + tricAlignment.getAAinSeq('J', resInThermoB) +
				" " + tricAlignment.getAAinSeq('J', resInThermoB,'I'));
		System.out.println("unitA " + tricAlignment.getAAinSeq('J', resInThermoB,'A') + 
				vec.getProfile('A', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitB " + tricAlignment.getAAinSeq('J', resInThermoB,'B') + 
				vec.getProfile('B', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitG " + tricAlignment.getAAinSeq('J', resInThermoB,'G') + 
				vec.getProfile('G', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitD " + tricAlignment.getAAinSeq('J', resInThermoB,'D') + 
				vec.getProfile('D', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitE " + tricAlignment.getAAinSeq('J', resInThermoB,'E') + 
				vec.getProfile('E', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitH " + tricAlignment.getAAinSeq('J', resInThermoB,'H') + 
				vec.getProfile('H', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitQ " + tricAlignment.getAAinSeq('J', resInThermoB,'Q') + 
				vec.getProfile('Q', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitZ " + tricAlignment.getAAinSeq('J', resInThermoB,'Z') + 
				vec.getProfile('Z', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println();	
		
		resInThermoB = 206;
		System.out.println("Therm " + resInThermoB + " " + tricAlignment.getAAinSeq('J', resInThermoB) +
				" " + tricAlignment.getAAinSeq('J', resInThermoB,'I'));
		System.out.println("unitA " + tricAlignment.getAAinSeq('J', resInThermoB,'A') + 
				vec.getProfile('A', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitB " + tricAlignment.getAAinSeq('J', resInThermoB,'B') + 
				vec.getProfile('B', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitG " + tricAlignment.getAAinSeq('J', resInThermoB,'G') + 
				vec.getProfile('G', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitD " + tricAlignment.getAAinSeq('J', resInThermoB,'D') + 
				vec.getProfile('D', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitE " + tricAlignment.getAAinSeq('J', resInThermoB,'E') + 
				vec.getProfile('E', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitH " + tricAlignment.getAAinSeq('J', resInThermoB,'H') + 
				vec.getProfile('H', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitQ " + tricAlignment.getAAinSeq('J', resInThermoB,'Q') + 
				vec.getProfile('Q', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println("unitZ " + tricAlignment.getAAinSeq('J', resInThermoB,'Z') + 
				vec.getProfile('Z', tricAlignment.getNewResNum('J', resInThermoB, 'A')-1));
		System.out.println();			
		
		int resInThermoA = 85;
		System.out.println("Therm " + resInThermoA + " " + tricAlignment.getAAinSeq('I', resInThermoA) +
				" " + tricAlignment.getAAinSeq('I', resInThermoA,'J'));
		System.out.println("unitA " + tricAlignment.getAAinSeq('I', resInThermoA,'A') + 
				vec.getProfile('A', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitB " + tricAlignment.getAAinSeq('I', resInThermoA,'B') + 
				vec.getProfile('B', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitG " + tricAlignment.getAAinSeq('I', resInThermoA,'G') + 
				vec.getProfile('G', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitD " + tricAlignment.getAAinSeq('I', resInThermoA,'D') + 
				vec.getProfile('D', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitE " + tricAlignment.getAAinSeq('I', resInThermoA,'E') + 
				vec.getProfile('E', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitH " + tricAlignment.getAAinSeq('I', resInThermoA,'H') + 
				vec.getProfile('H', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitQ " + tricAlignment.getAAinSeq('I', resInThermoA,'Q') + 
				vec.getProfile('Q', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitZ " + tricAlignment.getAAinSeq('I', resInThermoA,'Z') + 
				vec.getProfile('Z', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println();
		
		resInThermoA = 86;
		System.out.println("Therm " + resInThermoA + " " + tricAlignment.getAAinSeq('I', resInThermoA) +
				" " + tricAlignment.getAAinSeq('I', resInThermoA,'J'));
		System.out.println("unitA " + tricAlignment.getAAinSeq('I', resInThermoA,'A') + 
				vec.getProfile('A', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitB " + tricAlignment.getAAinSeq('I', resInThermoA,'B') + 
				vec.getProfile('B', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitG " + tricAlignment.getAAinSeq('I', resInThermoA,'G') + 
				vec.getProfile('G', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitD " + tricAlignment.getAAinSeq('I', resInThermoA,'D') + 
				vec.getProfile('D', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitE " + tricAlignment.getAAinSeq('I', resInThermoA,'E') + 
				vec.getProfile('E', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitH " + tricAlignment.getAAinSeq('I', resInThermoA,'H') + 
				vec.getProfile('H', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitQ " + tricAlignment.getAAinSeq('I', resInThermoA,'Q') + 
				vec.getProfile('Q', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitZ " + tricAlignment.getAAinSeq('I', resInThermoA,'Z') + 
				vec.getProfile('Z', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println();

		resInThermoA = 87;
		System.out.println("Therm " + resInThermoA + " " + tricAlignment.getAAinSeq('I', resInThermoA) +
				" " + tricAlignment.getAAinSeq('I', resInThermoA,'J'));
		System.out.println("unitA " + tricAlignment.getAAinSeq('I', resInThermoA,'A') + 
				vec.getProfile('A', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitB " + tricAlignment.getAAinSeq('I', resInThermoA,'B') + 
				vec.getProfile('B', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitG " + tricAlignment.getAAinSeq('I', resInThermoA,'G') + 
				vec.getProfile('G', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitD " + tricAlignment.getAAinSeq('I', resInThermoA,'D') + 
				vec.getProfile('D', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitE " + tricAlignment.getAAinSeq('I', resInThermoA,'E') + 
				vec.getProfile('E', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitH " + tricAlignment.getAAinSeq('I', resInThermoA,'H') + 
				vec.getProfile('H', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitQ " + tricAlignment.getAAinSeq('I', resInThermoA,'Q') + 
				vec.getProfile('Q', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitZ " + tricAlignment.getAAinSeq('I', resInThermoA,'Z') + 
				vec.getProfile('Z', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println();
		
		resInThermoA = 88;
		System.out.println("Therm " + resInThermoA + " " + tricAlignment.getAAinSeq('I', resInThermoA) +
				" " + tricAlignment.getAAinSeq('I', resInThermoA,'J'));
		System.out.println("unitA " + tricAlignment.getAAinSeq('I', resInThermoA,'A') + 
				vec.getProfile('A', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitB " + tricAlignment.getAAinSeq('I', resInThermoA,'B') + 
				vec.getProfile('B', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitG " + tricAlignment.getAAinSeq('I', resInThermoA,'G') + 
				vec.getProfile('G', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitD " + tricAlignment.getAAinSeq('I', resInThermoA,'D') + 
				vec.getProfile('D', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitE " + tricAlignment.getAAinSeq('I', resInThermoA,'E') + 
				vec.getProfile('E', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitH " + tricAlignment.getAAinSeq('I', resInThermoA,'H') + 
				vec.getProfile('H', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitQ " + tricAlignment.getAAinSeq('I', resInThermoA,'Q') + 
				vec.getProfile('Q', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitZ " + tricAlignment.getAAinSeq('I', resInThermoA,'Z') + 
				vec.getProfile('Z', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println();

		resInThermoA = 89;
		System.out.println("Therm " + resInThermoA + " " + tricAlignment.getAAinSeq('I', resInThermoA) +
				" " + tricAlignment.getAAinSeq('I', resInThermoA,'J'));
		System.out.println("unitA " + tricAlignment.getAAinSeq('I', resInThermoA,'A') + 
				vec.getProfile('A', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitB " + tricAlignment.getAAinSeq('I', resInThermoA,'B') + 
				vec.getProfile('B', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitG " + tricAlignment.getAAinSeq('I', resInThermoA,'G') + 
				vec.getProfile('G', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitD " + tricAlignment.getAAinSeq('I', resInThermoA,'D') + 
				vec.getProfile('D', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitE " + tricAlignment.getAAinSeq('I', resInThermoA,'E') + 
				vec.getProfile('E', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitH " + tricAlignment.getAAinSeq('I', resInThermoA,'H') + 
				vec.getProfile('H', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitQ " + tricAlignment.getAAinSeq('I', resInThermoA,'Q') + 
				vec.getProfile('Q', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitZ " + tricAlignment.getAAinSeq('I', resInThermoA,'Z') + 
				vec.getProfile('Z', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println();

		
		resInThermoA = 90;
		System.out.println("Therm " + resInThermoA + " " + tricAlignment.getAAinSeq('I', resInThermoA) +
				" " + tricAlignment.getAAinSeq('I', resInThermoA,'J'));
		System.out.println("unitA " + tricAlignment.getAAinSeq('I', resInThermoA,'A') + 
				vec.getProfile('A', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitB " + tricAlignment.getAAinSeq('I', resInThermoA,'B') + 
				vec.getProfile('B', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitG " + tricAlignment.getAAinSeq('I', resInThermoA,'G') + 
				vec.getProfile('G', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitD " + tricAlignment.getAAinSeq('I', resInThermoA,'D') + 
				vec.getProfile('D', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitE " + tricAlignment.getAAinSeq('I', resInThermoA,'E') + 
				vec.getProfile('E', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitH " + tricAlignment.getAAinSeq('I', resInThermoA,'H') + 
				vec.getProfile('H', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitQ " + tricAlignment.getAAinSeq('I', resInThermoA,'Q') + 
				vec.getProfile('Q', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitZ " + tricAlignment.getAAinSeq('I', resInThermoA,'Z') + 
				vec.getProfile('Z', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println();

		resInThermoA = 91;
		System.out.println("Therm " + resInThermoA + " " + tricAlignment.getAAinSeq('I', resInThermoA) +
				" " + tricAlignment.getAAinSeq('I', resInThermoA,'J'));
		System.out.println("unitA " + tricAlignment.getAAinSeq('I', resInThermoA,'A') + 
				vec.getProfile('A', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitB " + tricAlignment.getAAinSeq('I', resInThermoA,'B') + 
				vec.getProfile('B', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitG " + tricAlignment.getAAinSeq('I', resInThermoA,'G') + 
				vec.getProfile('G', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitD " + tricAlignment.getAAinSeq('I', resInThermoA,'D') + 
				vec.getProfile('D', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitE " + tricAlignment.getAAinSeq('I', resInThermoA,'E') + 
				vec.getProfile('E', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitH " + tricAlignment.getAAinSeq('I', resInThermoA,'H') + 
				vec.getProfile('H', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitQ " + tricAlignment.getAAinSeq('I', resInThermoA,'Q') + 
				vec.getProfile('Q', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitZ " + tricAlignment.getAAinSeq('I', resInThermoA,'Z') + 
				vec.getProfile('Z', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println();

		resInThermoA = 92;
		System.out.println("Therm " + resInThermoA + " " + tricAlignment.getAAinSeq('I', resInThermoA) +
				" " + tricAlignment.getAAinSeq('I', resInThermoA,'J'));
		System.out.println("unitA " + tricAlignment.getAAinSeq('I', resInThermoA,'A') + 
				vec.getProfile('A', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitB " + tricAlignment.getAAinSeq('I', resInThermoA,'B') + 
				vec.getProfile('B', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitG " + tricAlignment.getAAinSeq('I', resInThermoA,'G') + 
				vec.getProfile('G', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitD " + tricAlignment.getAAinSeq('I', resInThermoA,'D') + 
				vec.getProfile('D', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitE " + tricAlignment.getAAinSeq('I', resInThermoA,'E') + 
				vec.getProfile('E', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitH " + tricAlignment.getAAinSeq('I', resInThermoA,'H') + 
				vec.getProfile('H', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitQ " + tricAlignment.getAAinSeq('I', resInThermoA,'Q') + 
				vec.getProfile('Q', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println("unitZ " + tricAlignment.getAAinSeq('I', resInThermoA,'Z') + 
				vec.getProfile('Z', tricAlignment.getNewResNum('I', resInThermoA, 'A')-1));
		System.out.println();

	}

	
	
	
	
}
