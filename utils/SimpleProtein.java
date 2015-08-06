package utils;

import static utils.FileProcessor.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Created by Ziv_BA on 30/07/2015.
 */
public class SimpleProtein implements Iterable<SimpleProtein.ProtChain> {

    private File source;
    private List<ProtChain> protChains; // list of the seperate AminoAcid chains
    private List<String> hetAtmAndFooter;   // array of the remaining HeteroAtoms + footer tags.


    public File getSource() {
        return source;
    }

    /**
     * constructor for SimpleProtein from PDB file.</br>
     * creates atoms for each ATOM line, combine atoms by resSeq as residues, all residues of a specific
     * chain in a ProtChain object. each chain ends with a TER line. </br>
     * hetero atoms, master tag and end tag are all injected into the hetAtmAndFooter String array without
     * processing.
     * @param pdbFile File Object pointing to a PDB file.
     * @throws IOException if an error occurs while reading the file.
     */
    public SimpleProtein(File pdbFile) throws IOException {
        source = pdbFile;
        protChains = new ArrayList<>();
        // create a chain 2D array and the first chain list.
        List<ArrayList<String>> chains = new ArrayList<>(); chains.add(new ArrayList<String>());

        List<String> hetAtm = new ArrayList<>();
        int chainCounter = 0;

        List<String> PDBStrArr = Files.readAllLines(pdbFile.toPath(), Charset.defaultCharset());


        // iterate over the entire PDB file, populate the chain lists with ATOM lines, and hetAtom with HETATM lines.
        // for each TER entry, add a new list for that chain.
        for (int i = 0; i < PDBStrArr.size(); i++) {
            String lineInFile = PDBStrArr.get(i);
            if (lineInFile.startsWith("ATOM")) {
                chains.get(chainCounter).add(lineInFile);
            } else if (lineInFile.startsWith("TER")) {
                chains.get(chainCounter).add(lineInFile);
                chainCounter++;
                chains.add(new ArrayList<String>());
            }else if (lineInFile.matches("("+FOOTER_TAGS+").*")) {
                hetAtm.add(lineInFile);
            }

        }

        for (List<String> chain : chains) {
            if (chain != null && chain.size() >1) {
                protChains.add(new ProtChain(chain));
            }
        }

        hetAtmAndFooter = hetAtm;


    }

    /**
     * dump all the atom sequences and footer back to a PDB file.
     * @param destination
     * @throws IOException
     */
    public void writePDB(File destination) throws IOException {
        FileWriter FW = new FileWriter(destination);

        // write structural atoms
        for (ProtChain chain : protChains) {
            for (AminoAcid acid : chain){
                for (SimpleAtom atom : acid){
                    FW.write(atom.getOriginalString() + "\n");
                }
            }

        }
        //write HeteroAtoms and footer tags.
        for (String line : hetAtmAndFooter){
            FW.write(line);
        }
        FW.close();
    }

    public void createPermutations() throws IOException {

        File tempFolder = new File(source.getParent()+File.separator+"_temp");

        if (tempFolder.isDirectory()){
            System.out.println("temp folder already exists at:\n" + tempFolder.getAbsolutePath());
        } else {
            if (tempFolder.mkdir()) {
                System.out.println("Created temp processing folder \n"+ tempFolder.getAbsolutePath());
            } else {
                System.out.println("Temp Folder not created");
            }
        }

        ProteinActions.iterateAllAcids(this, tempFolder);
    }

    public void createPermutations(File tempFolder) throws IOException {
        if (tempFolder.isDirectory()){
            System.out.println("temp folder already exists at:\n" + tempFolder.getAbsolutePath());
        } else {
            if (tempFolder.mkdir()) {
                System.out.println("Created temp processing folder \n"+ tempFolder.getAbsolutePath());
            } else {
                System.out.println("Temp Folder not created");
            }
        }
        ProteinActions.iterateAllAcids(this, tempFolder);

    }


        @Override
    public Iterator<ProtChain> iterator() {
        return protChains.iterator();
    }


    /**
     * helper class to bulk residues together in respective chains.
     * also performs the actual processing from string array to molecular elements.
     */
    protected class ProtChain implements Iterable<AminoAcid>
    {
        private char chainID;
        private List<AminoAcid> residues = new ArrayList<>();

        /**
         * constructor creating a chain from a list of strings (assume all strings are for a single chain)
         * @param sourceList
         */
        protected ProtChain(List<String> sourceList)
        {
            List<String> tempAtomList = new ArrayList<>();
            int workingResSeq = Integer.valueOf(sourceList.get(0).substring(RES_SEQ_START,RES_SEQ_END+1)
                    .trim());

            for (int i = 0; i < sourceList.size(); i++) {

                String tempCurrLine = sourceList.get(i);
                int tempCurrResSeq = Integer.valueOf(tempCurrLine.substring(RES_SEQ_START, RES_SEQ_END + 1).trim());

                if (tempCurrResSeq == workingResSeq){
                    tempAtomList.add(tempCurrLine);
                } else {
                    workingResSeq = tempCurrResSeq;
                    residues.add(new AminoAcid(tempAtomList));
                    tempAtomList.clear();
                }

            }
            chainID = residues.get(0).getChainID();

        }

        public char getChainID() {
            return chainID;
        }

        protected void addRes(AminoAcid res) {
            residues.add(res);
        }

        @Override
        public Iterator<AminoAcid> iterator() {
            return residues.iterator();
        }
    }
}
