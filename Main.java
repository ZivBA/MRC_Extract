import utils.Scoring.MRC_Score;
import utils.fileUtilities.FileProcessor;
import utils.fileUtilities.MRC_Map_New;
import utils.molecularElements.ProteinActions;
import utils.molecularElements.SimpleProtein;

import java.io.File;
import java.io.IOException;


public class Main {

    /**
     * take command line argument 0 as path to MRC file, print out highest intensity + coords
     *
     * @param args path to MRC file
     */
    public static void main(String[] args) {
        String mrcpath = "//home//zivben//IdeaProjects//TestFilesForBBGenerator//3j06.mrc";
        FileProcessor FP = new FileProcessor(new File(args[0]), true);
        try {
            ProteinActions.stripAndAllALAToFile(FP.getSource(), FP.getDest());
            SimpleProtein processedProt = new SimpleProtein(FP.getDest());

            //            processedProt.createPermutations();
            MRC_Score scoreMap = new MRC_Score(new MRC_Map_New(mrcpath), processedProt);
            scoreMap.scoreProtein();
            scoreMap.calcZvalue();
            scoreMap.dispHist();

        } catch (IOException e) {
            e.printStackTrace();
        }

//


    }

}