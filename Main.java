import utils.FileProcessor;
import utils.ProteinActions;
import utils.SimpleProtein;

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
            ProteinActions.stripAndAllALA(FP.getSource(), FP.getDest());
            SimpleProtein processedProt = new SimpleProtein(FP.getDest());

            processedProt.createPermutations();


        } catch (IOException e) {
            e.printStackTrace();
        }

//


    }

}