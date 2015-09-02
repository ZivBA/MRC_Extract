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
	    String mrcpath = "//home//zivben//IdeaProjects//sputnikVirophage//EMD-5495.mrc";
	    //	    String mrcpath = "//home//zivben//IdeaProjects//TestFilesForBBGenerator//3j06.mrc";

	    SimpleProtein sourceProtein = null;
	    try {
		    if (args.length == 2) {
			    MRC_Score myScore = new MRC_Score(args[1], args[0]);
			    myScore.scoreProtein();
			    myScore.calcZvalue();
			    myScore.createCSVs();
		    } else if (args.length == 3) {
			    MRC_Score myScore = new MRC_Score(args[1], args[0], args[2]);
			    myScore.scoreProtein();
			    myScore.calcZvalue();
			    myScore.createCSVs();

		    }

	    } catch (IOException e) {
		    e.printStackTrace();
	    }


    }

}