package alignment;

import meshi.util.crossLinking.MySequence;
import meshi.util.crossLinking.MySequenceList;
import meshi.util.file.File2StringArray;

public class RvalAligner {
	
	public static double singleRvalAlignerRun(String seq , RvalSequence Rseq , boolean toPrintAlignment) {
		AminoAcidSequence AAseq = new AminoAcidSequence(seq);
//		SmithWatermanSolver SW = new SmithWatermanSolver(Rseq, AAseq, new RvalScoringScheme());
//		SW.printBestAlignment();
//		SW.mostMatchesAlignment(1000);
		NeedlemanWunchSolver NW = new NeedlemanWunchSolver(Rseq, AAseq, new RvalScoringScheme());
		if (toPrintAlignment) {
			NW.printAlignment();
			System.out.println();
		}
		return NW.alignmentScore();
	}
	
//
//	public static void main(String[] args) {
//		MySequenceList swissProt = new MySequenceList("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\uniprot-all_06_2015_allSeqs_in.fasta");
//		System.out.println("Total num of SwissProt seqs: " + swissProt.size());
//		String[] lines = File2StringArray.f2a("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\all_domains.txt");
//		for (int c=0 ; c<lines.length ; c++) {
//			MySequenceList mySeqList = new MySequenceList("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\PDBs\\"+lines[c]+".fasta");
//			RvalSequence Rseq = new RvalSequence("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\Run_not_minimized_with_HET_profiles\\"+lines[c]+"_Rval_profile.txt");
//			System.out.println("\n" + lines[c]  + ":  Structure positions: " + Rseq.size() + "   Length of *true* seq: " + mySeqList.get(0).seq().length() + "\n---------------------------------------------------------------------------------------");
//			double sameSeqScore = singleRvalAlignerRun(mySeqList.get(0).seq() , Rseq, true);
//			double largest = 0.0;
//			int indOflargest = -1;
//			int validSeqCounter = 0;
//			int worseScore = 0;
//			for (int ccc=0 ; ccc<swissProt.size() ; ccc++) {
//				MySequence mySeq = swissProt.get(ccc);
////				if (seqCounter%5000 == 0) {
////					System.out.println(seqCounter + ": " + worseScore + " with worse score.");
////				}
//				if (!(mySeq.seq().indexOf('X')>-1) &&
//						!(mySeq.seq().indexOf('B')>-1) &&
//						!(mySeq.seq().indexOf('J')>-1) &&
//						!(mySeq.seq().indexOf('Z')>-1) &&
//						!(mySeq.seq().indexOf('O')>-1) &&
//						!(mySeq.seq().indexOf('U')>-1)) {
//					if (mySeq.seq().length() >= Rseq.size()) {
//						double score = singleRvalAlignerRun(mySeq.seq() , Rseq , false);
//						if (score>-10000) {
//							validSeqCounter++;
//							if ((score>largest) & (score != sameSeqScore)) {
//								largest = score;
//								indOflargest = ccc;
//							}
//							if (score>sameSeqScore) {
//								worseScore++;
//							}
//						}
//					}
//				}
//			}
//			System.out.println(worseScore + " entries out of valid " + validSeqCounter + " are with better score than the *true*. The best score that is not the *true* sequence is: " + largest);
//			if (worseScore>0) {
//				MySequence mySeq = swissProt.get(indOflargest);
//				System.out.println("The best alignment is " + mySeq.title().substring(0,Math.min(mySeq.title().length()-1, 30)));
//				singleRvalAlignerRun(mySeq.seq() , Rseq, true);
//			}
//		}
//	}
//
//
//	public static void oldmain(String[] args) {
////		// Sequence of 4GWP_D
////		String S2 = "MSNQALYEKLEQTRTILSVKLAELINMTTIADRNDDDEGSFAQENSELAVATTSVMMVNNQTMQLIKNVQDLLILTRSIKEKWLLNQIPVTEHSKVTRFDEKQIEELLDNCIETFVAEKTT";
////		RvalSequence seq1 = new RvalSequence("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\For_Roger\\4GWP_D_Rval_profile.txt");
////		// Sequence of 4GWP_A
////		String S2 = "MQPPYIQERLKSLNDIETQLCSMLQEASQVTFIFGELKRGNESVKPQFENHVKQFYERLDKSTTQLRKEIQLLDENVGTRLLPINVNKKALGQDTEKMEEQLDLLSAILDPSKSK";
////		RvalSequence seq1 = new RvalSequence("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\For_Roger\\4GWP_A_Rval_profile.txt");
////		// Sequence of 4H54
////		String S2 = "AIKKTTEIDAILLNLNKAIDAHYQWLVSMFHSVVARDASKPEITDNHSYGLAQFGRWIDHLGPLDNDELPYVRLMDSAHQHMHNCGRELMLAIVENHWQDAHFDAFQEGLLSFTAALTDYKIYLLTIRSNMDVLTGLPGRRVLDESFDHQLRNAEPLNLYLMLLDIDRFKLVNDTYGHLIGDVVLRTLATYLASWTRDYETVYRYGGEEFIIIVKAANDEEACRAGVRICQLVDNHAITHSEGHINITVTAGVSRAFPEEPLDVVIGRADRAMYEGKQTGRNRCMFIDEQNVINRVLEHHHHHH";
////		RvalSequence seq1 = new RvalSequence("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\Run_not_minimized_with_HET_combined\\4H54_A_Rval_profile.txt");
////		// Sequence of 4LTP
////		String S2 = "GPSSPSLLRAIPGIAWIALLLLVIFYVFAVMGTKLFAQSFPEWFGTLGASMYTLFQVMTLESWSMGIARPVIEAYPWAWIYFVSFILVSSFTVLNLFIGIIIESMQSAHWEAEDAKRIEQEQRAHDERLEMLQLIRDLSSKVDRLERRSGKR";
////		RvalSequence seq1 = new RvalSequence("C:\\Users\\Nir\\Rvals\\4LTP_Rval_profile.txt");
////		// Sequence of 4IW0
////		String S2 = "GSQSTSNHLWLLSDILGQGATANVFRGRHKKTGDLFAIKVFNNISFLRPVDVQMREFEVLKKLNHKNIVKLFAIEEETTTRHKVLIMEFCPCGSLYTVLEEPSNAYGLPESEFLIVLRDVVGGMNHLRENGIVHRNIKPGNIMRVIGEDGQSVYKLTDFGAARELEDDEQFVSLYGTEEYLHPDMYERAVLRKDHQKKYGATVDLWSIGVTFYHAATGSLPFRPFEGPRRNKEVMYKIITGKPSGAISGVQKAENGPIDWSGDMPVSCSLSRGLQVLLTPVLANILEADQEKCWGFDQFFAETSDILHRMVIHVFSLQQMTAHKIYIHSYNTATIFHELVYKQTKIISSNQELIYEGRRLVLEPGRLAQHFPKTTEENPIFVVSREPLNTIGLIYEKISLPKVHPRYDLDGDASMAKAITGVVCYACRIASTLLLYQELMRKGIRWLIELIKDDYNETVHKKTEVVITLDFCIRNIEKTVKVYEKLMKINLEAAELGEISDIHTKLLRLSSSQGTIETSLQDIDSRLSPGGSLADAWAHQEGTHPKDRNVEKLQVLLNCMTEIYYQFKKDKAERRLAYNEEQIHKFDKQKLYYHATKAMTHFTDECVKKYEAFLNKSEEWIRKMLHLRKQLLSLTNQCFDIEEEVSKYQEYTNELQET";
////		RvalSequence seq1 = new RvalSequence("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\Run_not_minimized_with_HET_combined\\4IW0_A_Rval_profile.txt");
////		// Sequence of 4HYD
////		String S2 = "MQIRDWLPLLGMPLMLLFVQIIAIVLVMPMQAAGLVAFENPSSVANPLIFIGMLLAFTLVLLVLLRTGGRRFIAAFIGFALFMTFLYIFGALSLLALGPTTAAAAGTLIGAVAVTALLYLYPEWYVIDILGVLISAGVASIFGISLEPLPVLVLLVLLAVYDAISVYRTKHMITLAEGVLETKAPIMVVVPKRADYSFRKEGLNISEGEERGAFVMGMGDLIMPSILVVSSHVFVDAPAVLWTLSAPTLGAMVGSLVGLAVLLYFVNKGNPQAGLPPLNGGAILGFLVGAALAGSFSWLPF";
////		RvalSequence seq1 = new RvalSequence("C:\\Users\\Nir\\Rvals\\4HYD_Rval_profile.txt");
//
////		// Sequence of 3QIL_B
////		String S2 = "MGSSHHHHHHSSGLVPRGSHMWKQSVELAKKDSLYKDAMQYASESKDTELAEELLQWFLQEEKRECFGACLFTCYDLLRPDVVLELAWRHNIMDFAMPYFIQVMKEYLTKVDKLDASESLRKEEE";
////		RvalSequence seq1 = new RvalSequence("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\TEST\\3QIL_B_Rval_profile.txt");
//
////		// Sequence of 3ZZI_A
////		String S2 = "MGHHHHHHVSSTNGFSATRSTVIQLLNNISTKREVEQYLKYFTSVSQQQFAVIKVGGAIISDNLHELASCLAFLYHVGLYPIVLHGTGPQVNGRLEAQGIEPDYIDGIRITDEHTMAVVRKCFLEQNLKLVTALEQLGVRARPITSGVFTADYLDKDKYKLVGNIKSVTKEPIEASIKAGALPILTSLAETASGQMLNVNADVAAGELARVFEPLKIVYLNEKGGIINGSTGEKISMINLDEEYDDLMKQSWVKYGTKLKIREIKELLDYLPRSSSVAIINVQDLQKELFTDSGAGTMIRRGYKLVKRSSIGEFPSADALRKALQRDAGISSGKESVASYLRYLENSDFVSYADEPLEAVAIVKKDTNVPTLDKFVCSDAAWLNNVTDNVFNVLRRDFPALQWVVSENDANIAWHFDKSQGSYLKGGKVLFWYGIDDINTISELVENFVKSCDTASTLNSSASS";
////		RvalSequence seq1 = new RvalSequence("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\TEST\\3ZZI_A_Rval_profile.txt");
//
////		// Sequence of 3N7N_E
////		String S2 = "MTTLLQLLSNYYKAKLDSERIYNEYVQSQYEFASLDKPKKVVDETLFLQRQIAQLNKQLQLSFQENEKLLSVQKNQKALY";
////		RvalSequence seq1 = new RvalSequence("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\TEST\\3N7N_E_Rval_profile.txt");
//
//		// Sequence of 3N7N_A
//		String S2 = "MDPLTVYKNSVKQQIDSADLLVANLVNENFVLSEKLDTKATEIKQLQKQIDSLNAQVKELKTQTSQQAENSEVIKDLYEYLCNVRVHKSYEDDSGLWFDISQGTHSGGSSDDYSIMDYKLGFVKGQAQVTEVIYAPVLKQRSTEELYSLQSKLPEYLFETLSFPLSSLNQFYNKIAKSLNKKREKKDETE";
//		RvalSequence seq1 = new RvalSequence("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\TEST\\3N7N_A_Rval_profile.txt");
//
//
////		// Sequence of 3W3A_A
////		String S2 = "MIQGVIQKIAGPAVIAKGMLGARMYDICKVGEEGLVGEIIRLDGDTAFVQVYEDTSGLKVGEPVVSTGLPLAVELGPGMLNGIYDGIQRPLERIREKTGIYITRGVVVHALDREKKWAWTPMVKPGDEVRGGMVLGTVPEFGFTHKILVPPDVRGRVKEVKPAGEYTVEEPVVVLEDGTELKMYHTWPVRRARPVQRKLDPNTPFLTGMRILDVLFPVAMGGTAAIPGPFGSGKTVTQQSLAKWSNADVVVYVGCGERGNEMTDVLVEFPELTDPKTGGPLMHRTVLIANTSNMPVAAREASIYVGVTIAEYFRDQGFSVALMADSTSRWAEALREISSRLEEMPAEEGYPPYLAARLAAFYERAGKVITLGGEEGAVTIVGAVSPPGGDMSEPVTQSTLRIVGAFWRLDASLAFRRHFPAINWNGSYSLFTSALDPWYRENVAEDYPELRDAISELLQREAGLQEIVQLVGPDALQDAERLVIEVGRIIREDFLQQNAYHEVDAYCSMKKAYGIMKMILAFYKEAEAAIKRGVSIDEILQLPVLERIGRARYVSEEEFPAYFEEAMKEIQGAFKAL";
////		RvalSequence seq1 = new RvalSequence("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\TEST\\3W3A_A_Rval_profile.txt");
//
////		// Sequence of 4HULA_A
////		String S2 = "LDRFSFSVFLKEIRLLTALALPMLLAQVAQVGIGFVDTVMAGGAGKEDLAAVALGSSAFATVYITFMGIMAALNPMIAQLYGAGKTGEAGETGRQGIWFGLILGIFGMILMWAAITPFRNWLTLSDYVEGTMAQYMLFTSLAMPAAMVHRALHAYASSLNRPRLIMLVSFAAFVLNVPLNYIFVYGKFGMPALGGAGCGVATMAVFWFSALALWIYIAKEKFFRPFGLTAKFGKPDWAVFKQIWKIGAPIGLSYFLEASAFSFIVFLIAPFGEDYVAAQQVGISLSGILYMIPQSVGSAGTVRIGFSLGRREFSRARYISGVSLVSGWVLAVITVLSLVLFRSPLASMYNDDPAVLSIASTVLLFAGLFQPADFTQCIASYALRGYKVTKVPMFIHAAAFWGCGLLPGYLLAYRFDMGIYGFWTALIASLTIAAVALVWCLEKYSMELVKSHKAVSSGL";
////		RvalSequence seq1 = new RvalSequence("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\TEST\\4HUL_A_Rval_profile.txt");
//
////		// Sequence of 4HULA_A
////		String S2 = "GSPEFVQGNTDEAQEELAWKIAKMIVSDVMQQAQYDQPLEKSTKL";
////		RvalSequence seq1 = new RvalSequence("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\TEST\\3TMH_D_Rval_profile.txt");
//
////		// Sequence of 3T1P_A
////		String S2 = "NKITPNLAEFAFSLYRQLAHQSNSTNIFFSPVSIATAFAMLSLGTKADTHDEILEGLNFNLTEIPEAQIHEGFQELLRTLNQPDSQLQLTTGNGLFLSEGLKLVDKFLEDVKKLYHSEAFTVNFGDTEEAKKQINDYVEKGTQGKIVDLVKELDRDTVFALVNYIFFKGKWERPFEVKDTEEEDFHVDQVTTVKVPMMKRLGMFNIQHSKKLSSWVLLMKYLGNATAIFFLPDEGKLQHLENELTHDIITKFLENEDRRSASLHLPKLCITGTYDLKSVLGQLGITKVFSNGADLSGVTEEAPLKLSKAVHKAVLCIDEKGTEAAGAMFLEAIPRSIPPEVKFNKPFVFLMIEQNTKSPLFMGKVVNPTQK";
////		RvalSequence seq1 = new RvalSequence("C:\\Users\\Nir\\Check_R_Val\\Large_Scale_Alignment\\TEST\\3T1P_A_Rval_profile.txt");
//
//
//		AminoAcidSequence seq2 = new AminoAcidSequence(S2);
//		SmithWatermanSolver SW = new SmithWatermanSolver(seq1, seq2, new RvalScoringScheme());
//		double bestScore = SW.bestAlignmentScore();
//		SW.printBestAlignment();
//
//		MySequenceList seqList = new MySequenceList("C:\\Users\\Nir\\Check_R_Val\\Scan_a_single_residue\\uniprot_sprot.fasta");
//		System.out.println("--Number of sequences for threading: " + seqList.size());
//		double largest = 0.0;
//		int seqCounter = 0;
//		int worseScore = 0;
//		for (MySequence mySeq : seqList) {
//			seqCounter++;
//			if (seqCounter%50000 == 0) {
//				System.out.println(seqCounter + ": " + worseScore + " with worse score.");
//			}
//			if (!(mySeq.seq().indexOf('X')>-1) &&
//					!(mySeq.seq().indexOf('B')>-1) &&
//					!(mySeq.seq().indexOf('J')>-1) &&
//					!(mySeq.seq().indexOf('Z')>-1) &&
//					!(mySeq.seq().indexOf('O')>-1) &&
//					!(mySeq.seq().indexOf('U')>-1)) {
//				seq2 = new AminoAcidSequence(mySeq.seq());
//			 	SW = new SmithWatermanSolver(seq1, seq2, new RvalScoringScheme());
//				double score = SW.bestAlignmentScore();
//				if (score>largest) {
//					largest = score;
//				}
//				if (score>bestScore) {
//					System.out.println(mySeq.title());
//					SW.printBestAlignment();
//					worseScore++;
//				}
//			}
//		}
//		System.out.println(largest);
//		System.out.println(worseScore + " with worse score.");
//	}
	
	
	
	/**
	 * The code I used to add the sequences that were not in swissprot
 		String addTo = "";
		String res = lines[c]+" match not found!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
		boolean found = false;
		for (int sp=0 ; sp<swissProt.size() ; sp++) {
			if (mySeqList.get(0).seq().equals(swissProt.get(sp).seq())) {
				res = lines[c]+" match found: " + swissProt.get(sp).title();
				found = true;
			}
		}
		System.out.println(res);
		if (!found) {
			addTo += (">" +lines[c]+"-I added 19/06/15 \n" + mySeqList.get(0).seq() + "\n");
		}
		System.out.println("\n\n" + addTo);

	 **/

}
