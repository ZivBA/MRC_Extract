package meshi.molecularElements.allAtoms;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.GridCell;
import meshi.geometry.MatrixRow;
import meshi.molecularElements.Atom;
import meshi.molecularElements.AtomList;
import meshi.util.UpdateableException;

public class DistanceMatrixWhichUpdateOnlyMovingAtoms extends DistanceMatrix {
    private double[][] previousCoordinates;
    private int numberOfUpdates = 0;
    public DistanceMatrixWhichUpdateOnlyMovingAtoms(AtomList atomList, double rMax, int bondedListDepth) {
	super(atomList,rMax,0.0,bondedListDepth);
	initPreviousCoordinates();
    }

    public void update(int numberOfUpdates) throws UpdateableException {
	if (numberOfUpdates == this.numberOfUpdates+1) {
	    this.numberOfUpdates++;
	    updateMovingAtoms();
            
	}
	else if (numberOfUpdates != this.numberOfUpdates) 
	    throw new RuntimeException("Something weird with DistanceMatrix.update(int numberOfUpdates)\n"+
				       "numberOfUpdates = "+numberOfUpdates+" this.numberOfUpdates = "+this.numberOfUpdates);
    }				       

    private void updateMovingAtoms()  throws UpdateableException{
	int length = matrix.length;
	for (int iRow = 0; iRow <length; iRow++)
	    if (moved(matrix[iRow].atom,iRow)) {
		matrix[iRow].updateIncludingMirrors();
	    }
	int size = atomList.size();
        grid.build();
        for (int iatom = 0; iatom < size; iatom++) {
	    MatrixRow row = matrix[iatom];
	    GridCell gridCell = grid.getCell(row.atom);
	    row.addCell(gridCell);
        }
	getNonBondedList();
    }
    private boolean moved(Atom atom, int index) {
	double[] prev = previousCoordinates[index];
	if ((atom.x() != prev[0]) ||
	    (atom.y() != prev[1]) ||
	    (atom.z() != prev[2])) {
	    prev[0] = atom.x();
	    prev[1] = atom.y();
	    prev[2] = atom.z();
	    return true;
	}
	return false;
    }

    public void initPreviousCoordinates() {
	int size = atomArray.length;
	previousCoordinates = new double[size][3];
	for (int i = 0; i< size; i++) {
	    Atom atom = (Atom) atomArray[i];
	    previousCoordinates[i][0] = atom.x();
	    previousCoordinates[i][1] = atom.y();
	    previousCoordinates[i][2] = atom.z();
	}
    }
}
