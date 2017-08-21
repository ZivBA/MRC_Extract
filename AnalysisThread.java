
import utils.ScoreUtilities.ScoringGeneralHelpers;
import utils.Scoring.MRC_Score;
import utils.molecularElements.AminoAcid;
import utils.molecularElements.SimpleAtom;
import utils.molecularElements.SimpleProtein;

import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Callable;

/**
 * Created by zivben on 11/12/16.
 */
public class AnalysisThread implements Callable<HashMap<Character,HashMap<String,List<Double[]>>>> {
	
	private final String[] args;
	
	public AnalysisThread(String[] tempArgs) {
		this.args = tempArgs;
	}
	
	@Override
	public HashMap<Character,HashMap<String,List<Double[]>>> call() throws Exception{
		
		HashMap<Character,HashMap<String,List<Double[]>>> resultMap = new HashMap<>();
		for (char aAcid : ScoringGeneralHelpers.singleLetters){
			resultMap.put(aAcid, new HashMap<String,List<Double[]>>());
		}
		
		MRC_Score myScore = null;
		try {
			// if argument string does not include chain to process, then entire protein needs to be processed (not used).
			if (args.length == 2) {
				myScore = new MRC_Score(args[1], args[0]);
				char[] chains = myScore.getMyProt().getChains();
				
				for (char chain : chains) {
					myScore = new MRC_Score(myScore.getMyProt(), String.valueOf(chain), myScore.getMyMap());
					SimpleProtein.ProtChain workingChaing = myScore.getMyProt().getChain(chain);
					for (AminoAcid acid : workingChaing){
						float[] alphaCoords = acid.getCalpha().getAtomCoords();
						double cAlphaScore = myScore.getMyMap().val(alphaCoords[0],alphaCoords[1],alphaCoords[2]);
						for (SimpleAtom atom : acid){
							float[] coords = atom.getAtomCoords();
							
							try {
								atom.setAtomScore(myScore.getMyMap().val(coords[0], coords[1], coords[2]));
								
							} catch (RuntimeException e) {
								atom.setAtomScore(0.0);
							}
							if (!resultMap.get(acid.getSingleLetter()).containsKey(atom.getName())){
								resultMap.get(acid.getSingleLetter()).put(atom.getName(),new LinkedList<Double[]>());
								resultMap.get(acid.getSingleLetter()).get(atom.getName()).add(new Double[]{cAlphaScore, atom.getAtomScore()});
							} else{
								resultMap.get(acid.getSingleLetter()).get(atom.getName()).add(new Double[]{cAlphaScore, atom.getAtomScore()});
							}
						}
					}
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return resultMap;
		
		
	}
	
}

