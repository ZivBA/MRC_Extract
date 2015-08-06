package meshi.sequences;
import meshi.molecularElements.Residue;
/**
 * A pair of residues each from a different Protein object and a string indicating which 
 * atoms in these residues should be aligned. 
 * Represents a column in a sequence alignment (currently no multiple sequences).
 **/ 

public class ResidueAlignmentColumn extends AlignmentColumn{
    protected ResidueAlignmentColumn(int numberOfRows){
	super(numberOfRows);
    }

    public ResidueAlignmentColumn(Residue residue0, Residue residue1) {
	super(new ResidueAlignmentCell(residue0),
	      new ResidueAlignmentCell(residue1));
    }
    public ResidueAlignmentColumn(ResidueAlignmentCell cell0, ResidueAlignmentCell cell1) {
	super(cell0, cell1);
    }
    
    public String toString() {
	return ""+residue0()+"\t"+residue1();
    }
    
    public Residue residue0() { return (Residue) cell0().obj;}
    public Residue residue1() { return (Residue) cell1().obj;}
}
