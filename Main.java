import utils.Scoring.MRC_Score;
import utils.molecularElements.SimpleProtein;

import java.io.IOException;


public class Main {

    /**
     * take command line argument 0 as path to MRC file, print out highest intensity + coords
     *
     * @param args path to MRC file
     */
    public static void main(String[] args) {
        String mrcpath = "//home//zivben//IdeaProjects//TestFilesForBBGenerator//3j06.mrc";

	    SimpleProtein sourceProtein = null;
	    try {
		    MRC_Score myScore = new MRC_Score(mrcpath, args[0]);
		    myScore.scoreProtein();
		    myScore.calcZvalue();
		    System.out.println();


	    } catch (IOException e) {
		    e.printStackTrace();
	    }


    }

}