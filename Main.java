import utils.Scoring.MRC_Score;
import utils.fileUtilities.FileProcessor;

import java.io.File;


public class Main {

    /**
     * take command line argument 0 as path to MRC file, print out highest intensity + coords
     *
     * @param args path to MRC file
     */
    public static void main(String[] args) {
        String mrcpath = "//home//zivben//IdeaProjects//TestFilesForBBGenerator//3j06.mrc";


        FileProcessor FP = new FileProcessor(new File(args[0]), true);
        MRC_Score scoreMap = MRC_Score.StartFromScratch(FP, mrcpath);

        /*try {
            ProteinActions.stripAndAllALAToFile(FP.getSource(), FP.getDest());
            SimpleProtein processedProt = new SimpleProtein(FP.getDest());


            MRC_Score scoreMap = new MRC_Score(new MRC_Map_New(mrcpath), processedProt);
            scoreMap.scoreProtein();
            scoreMap.calcZvalue();
            scoreMap.dispHist();

        } catch (IOException e) {
            e.printStackTrace();
        }*/

//


    }

}