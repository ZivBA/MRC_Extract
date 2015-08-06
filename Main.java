import utils.fileUtilities.FileProcessor;
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

        FileProcessor FP = new FileProcessor(new File(args[0]), true);
        try {
            ProteinActions.stripAndAllALAToFile(FP.getSource(), FP.getDest());
            SimpleProtein processedProt = new SimpleProtein(FP.getDest());

            processedProt.createPermutations();


        } catch (IOException e) {
            e.printStackTrace();
        }

//


    }

}