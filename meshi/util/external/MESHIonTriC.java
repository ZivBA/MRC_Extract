package meshi.util.external;

import meshi.geometry.DistanceMatrix;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.rotamericTools.RotamericTools;

public class MESHIonTriC {
	
	private Protein complex = null;
	private DunbrackLib lib = null;
	
	public MESHIonTriC(String complexPDB, DunbrackLib lib) {
		complex = ComplexMESHIconversion.complex2meshi(complexPDB);
		this.lib = lib;
	}

	public double[][] putIntoNearestRot() {
		return RotamericTools.putIntoNearestRot(complex, new DistanceMatrix(complex.atoms(), 2.0, 1.0), lib,true);	
	}

	public double[][] getNearestRotInfo() {
		return RotamericTools.putIntoNearestRot(complex, new DistanceMatrix(complex.atoms(), 2.0, 1.0), lib,false);	
	}

	public void writeComplexToDisk(String fileName) {
		ComplexMESHIconversion.writeMEHSI2complex(complex, fileName);
	}
	

	public static void main(String[] args) {
		MeshiProgram.initRandom(0);
		CommandList commands = new CommandList(args[0]);
		MESHIonTriC meshiOnTriC  = new MESHIonTriC(args[1], new DunbrackLib(commands,0.99,100));	
		meshiOnTriC.putIntoNearestRot();
		meshiOnTriC.writeComplexToDisk(args[2]);
	}

}
